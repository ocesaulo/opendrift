  # This file is part of OpenDrift.
#
# OpenDrift is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2
#
# OpenDrift is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with OpenDrift.  If not, see <https://www.gnu.org/licenses/>.
#
# Copyright 2020, Knut-Frode Dagestad, MET Norway

"""
SedimentDrift is an OpenDrift module for drift and settling of sediments.
Based on work by Simon Weppe, MetOcean Solutions Ltd.
"""

import numpy as np
import logging; logger = logging.getLogger(__name__)
from opendrift.models.oceandrift import OceanDrift
from opendrift.models.oceandrift import Lagrangian3DArray
from opendrift.config import CONFIG_LEVEL_ESSENTIAL, CONFIG_LEVEL_BASIC, CONFIG_LEVEL_ADVANCED
from datetime import datetime
from scipy.integrate import solve_ivp
from scipy import sparse
from opendrift.models.bblm_sg2000 import sg2000_solve
from opendrift.models import floc_strength
from opendrift.models import rouse
from opendrift.models import settling
import math


def particle_motion_vectorized(t, y_flat, params):
    n = y_flat.size // 2
    y = y_flat.reshape((2, n))
    z = y[0, :]
    w = y[1, :]

    Cl, rho_f, rho_s, u_star, r, mu, g = params
    z_fixed = np.where(z == 0, 1e-8, z)

    lift_term = (0.5 * 3 * Cl * rho_f / (4 * r * rho_s)) * (u_star**2 * (np.log(2))**2) / (z_fixed**2)
    gravity_term = g * (1 - rho_f / rho_s)
    drag_term = (9 * mu / (2 * r**2 * rho_s)) * w

    dw_dt = lift_term - gravity_term - drag_term
    dz_dt = w
    dydt = np.vstack((dz_dt, dw_dt))
    return dydt.flatten()


def _resuspension_ode_jacobian_pattern(n):
    """Constant (row, col) index arrays for the block-structured Jacobian of the
    resuspension ODE (state ordered [z_0..z_{n-1}, w_0..w_{n-1}]).

    Non-zeros: d(dz/dt)/dw = I, d(dw/dt)/dz = diag, d(dw/dt)/dw = diag -> 3n entries.
    The system is decoupled across particles, so the Jacobian is sparse and the
    sparse BDF solve scales ~O(n) instead of the O((2n)^3) of a dense Jacobian.
    """
    idx = np.arange(n)
    rows = np.concatenate([idx, n + idx, n + idx])
    cols = np.concatenate([n + idx, idx, n + idx])
    return rows, cols


class SedimentElement(Lagrangian3DArray):
    variables = Lagrangian3DArray.add_variables([
        ('settled', {'dtype': np.uint8,  # 0 is active, 1 is settled
                     'units': '1',
                     'default': 0}),
        ('terminal_velocity', {'dtype': np.float32,
                               'units': 'm/s',
                               'default': -0.001}),  # 1 mm/s negative buoyancy
        # use the terminal_velocity_default to restore downwards terminal velocity after
        # particle has been given an upwards one for resuspension
        # ('terminal_velocity_default', {'dtype': np.float32,
        #                        'units': 'm/s',
        #                        'default': -0.001}),  
        ('grain_diameter', {'dtype': np.float32,
                 'units': 'm',
                 'default': 4e-6}),
        # used to give particle upward resuspension velocity for 1 timestep
        # ('counter', {'dtype': np.uint8,
        #               'units': '1',
        #               'default': 0}),
        # Distance from particle to nearest cell center above, value set to 99999
        # if particle not settled
        # ('dz_bot', {'dtype': np.float32,
        #               'units': 'm',
        #               'default': 99999.0}), 
        ('C_l', {'dtype': np.float32,
                      'units': 'unitless',
                      'default': 0.5}),
        ('tau_crit', {'dtype': np.float32,
                      'units': 'Pa',
                      'default': 0.09}),
        # Sediment behaviour class for the dynamic tau_crit closure:
        # 0 = non-cohesive/solid grain (Shields), 1 = cohesive floc (fractal strength).
        # Default 0; a size-cutoff default is applied at the first resuspension
        # step for elements left at 0 with a sub-cutoff grain size (see
        # _critical_shear_stress). Only used when tau_crit_mode != 'constant'.
        ('sed_class', {'dtype': np.uint8,
                      'units': '1',
                      'default': 0}),
        # Floc fractal dimension Df (cohesive elements). Marine flocs ~1.7-2.3.
        ('fractal_dim', {'dtype': np.float32,
                      'units': '1',
                      'default': 2.0}),
        # Fresh (as-deposited) solids volume fraction phi0 of a cohesive floc;
        # consolidation grows it toward consolidated_solids_fraction over time.
        ('phi0', {'dtype': np.float32,
                      'units': '1',
                      'default': 0.1}),
        # Time (s) an element has been settled on the bed; drives consolidation
        # of the cohesive critical shear stress. Reset to 0 on (re)settling and
        # on resuspension.
        ('time_since_settled', {'dtype': np.float32,
                      'units': 's',
                      'default': 0}),
        ('rho_s', {'dtype': np.float32,
                      'units': 'kgm-3',
                      'default': 2000}),
        # ('porosity', {'dtype': np.float32,
        #               'units': 'unitless',
        #               'default': 0.9}),
        # Molecular viscosity used for Stokes Law terminal velocity calculation.
        # Since this is a property of seawater, it would make sense for this
        # parameter to be stored in the model config. However, by putting it
        # in this location, the move_elements function has access to it, and
        # the move_elements function is where terminal velocity is calculated.
        ('viscosity_molecular', {'dtype': np.float32,
                      'units': 'kgm-1s-1',
                      'default': 1.4e-3}),
        # whether or not to use Stokes Law calculation or empirical value for terminal velocity
        ('use_stokes', {'dtype': np.uint8,
                      'units': '1',
                      'default': 1}),
        # keeps track of how many times a particle was resuspended
        ('times_resuspended', {'dtype': np.uint8,
                      'units': '1',
                      'default': 0}),
        ('beached', {'dtype': np.uint8,  # 0 is active, 1 is settled
                     'units': '1',
                     'default': 0}),
        ('latest_resuspension_height', {'dtype': np.float32,
                      'units': 'm',
                      'default': 0}),
        # ----- settling-velocity (terminal velocity) closure inputs -----
        # Primary (constituent) particle size of a fractal aggregate, used by the
        # 'maggi'/'maggi_permeable' settling models for the mass-size relation
        # (d/d0)**(Df-3). For a solid grain set d0 = grain_diameter (the seeding
        # layer does this), giving an effective density == rho_s.
        ('d0', {'dtype': np.float32,
                      'units': 'm',
                      'default': 4e-6}),
        # Corey shape factor S/sqrt(L*I) of the grain (1 = sphere), used by the
        # 'bb16' settling model's drag-correction factors.
        ('corey_shape_factor', {'dtype': np.float32,
                      'units': '1',
                      'default': 1.0}),
        # Set to 1 once update_terminal_velocity() has computed this element's
        # settling velocity, so the static (non-dynamic) closure computes it once
        # at first activation instead of every step. Ignored when settling_dynamic.
        ('vel_set', {'dtype': np.uint8,
                      'units': '1',
                      'default': 0})
        ])

    def move_elements(self, other, indices):
        super(Lagrangian3DArray, self).move_elements(other, indices)
        # Set terminal velocity to something calculated using Stokes Law
        grain_diameter = other.grain_diameter
        # This viscosity is in [kg/(ms)] so do not need to multiply by density seawater
        viscosity = other.viscosity_molecular
        gravity = -9.81
        # Since this method is a part of the Lagrangian3DArray class, it does not have access
        # to the sea_water_density function from physics methods, thus I just chose a default value.
        rho_ocean = 1026.95
        rho_sed = other.rho_s
        term_vel = 2 / 9 * (rho_sed - rho_ocean) * gravity / viscosity * (grain_diameter / 2)**2 
        not_stokes = 1 - other.use_stokes
        term_vel = term_vel * other.use_stokes + other.terminal_velocity * not_stokes
        other.terminal_velocity = term_vel
        # other.terminal_velocity_default = term_vel
      


class SedimentDrift(OceanDrift):
    """Model for sediment drift, under development
    """

    ElementType = SedimentElement

    required_variables = {
        'x_sea_water_velocity': {'fallback': 0},
        'y_sea_water_velocity': {'fallback': 0},
        'sea_surface_height': {'fallback': 0},
        'upward_sea_water_velocity': {'fallback': 0},
        'x_wind': {'fallback': 0},
        'y_wind': {'fallback': 0},
        'sea_surface_wave_stokes_drift_x_velocity': {'fallback': 0},
        'sea_surface_wave_stokes_drift_y_velocity': {'fallback': 0},
        'sea_surface_wave_significant_height': {'fallback': 0},
        'sea_surface_wave_period_at_variance_spectral_density_maximum': {'fallback': 0},
        'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment': {'fallback': 0},
        'land_binary_mask': {'fallback': None},
        'sea_water_temperature': {'fallback': 10},
        'sea_water_salinity': {'fallback': 35},
        'ocean_vertical_diffusivity': {'fallback': 0.02,
                                      'profiles': True},
        'ocean_mixed_layer_thickness': {'fallback': 50},
        'sea_floor_depth_below_sea_level': {'fallback': 10000},
        }

    def __init__(self, *args, **kwargs):
        """ Constructor of SedimentDrift module
        """

        super(SedimentDrift, self).__init__(*args, **kwargs)

        self._add_config({
            'vertical_mixing:resuspension_threshold': {
                'type': 'float',
                'default': 0.2,
                'min': 0,
                'max': 3,
                'units': 'm/s',
                'description':
                'Sedimented particles will be resuspended if bottom current shear exceeds this value.',
                'level': CONFIG_LEVEL_ESSENTIAL
            }})

        # Currently this config value is not being used as it is also saved in the
        # SedimentELement object. However, if the calculation for terminal velocity
        # were to be moved out of the SedimentElement class, it would be useful to
        # store molecular viscosity here.
        self._add_config({
            'environment:molecular_viscosity': {
                'type': 'float',
                'default': 1.4e-3,
                'min': 1e-3,
                'max': 2e-3,
                'units': 'kg(ms)-1',
                'description':
                'Sedimented particles will be resuspended if bottom current shear exceeds this value.',
                'level': CONFIG_LEVEL_ESSENTIAL
            }})

        # ----- Settling (terminal) velocity closure -----
        # How each element's terminal_velocity is obtained (see models/settling.py):
        #   'prescribed' - LEGACY default: terminal_velocity is left as seeded
        #                  (the use_stokes/empirical path in SedimentElement.move_elements);
        #                  update_terminal_velocity() is a no-op. Reproduces old runs.
        #   'stokes'     - model-side Stokes law using the real local fluid density.
        #   'dietrich'   - Dietrich (1982) drag (sub-spherical grain).
        #   'bb16'       - Bagheri & Bonadonna (2016), uses corey_shape_factor.
        #   'maggi' / 'maggi_permeable' - Maggi (2013) fractal aggregate, uses d0,
        #                  fractal_dim.
        self._add_config({
            'vertical_mixing:settling_model': {
                'type': 'enum',
                'enum': ['prescribed', 'stokes', 'dietrich', 'bb16',
                         'maggi', 'maggi_permeable'],
                'default': 'prescribed',
                'description':
                'Closure for the per-element settling (terminal) velocity. '
                "Default 'prescribed' keeps the legacy seeded value; the other "
                'modes compute it from grain size/density/shape and the ambient '
                'fluid via models/settling.py.',
                'level': CONFIG_LEVEL_BASIC
            }})
        # When False (default) a non-'prescribed' settling velocity is computed once,
        # at each element's first active step, and then held fixed (cheap, seed-time
        # static). When True it is recomputed every step from the local T/S-derived
        # fluid density and viscosity -> dynamic terminal velocity, no other change.
        self._add_config({
            'vertical_mixing:settling_dynamic': {
                'type': 'bool',
                'default': False,
                'description':
                'Recompute the settling velocity every step from the local fluid '
                'properties (dynamic) instead of once at first activation (static).',
                'level': CONFIG_LEVEL_BASIC
            }})


        # Bottom boundary layer scheme used to compute the bed shear stress that
        # drives resuspension. 'sg2000' is the Styles & Glenn (2000) combined
        # wave-current model (reduces to a current-only log-law where no waves);
        # 'legacy' is the old c_d drag + diffusivity estimate.
        self._add_config({
            'vertical_mixing:bbl_scheme': {
                'type': 'enum',
                'enum': ['sg2000', 'legacy'],
                'default': 'sg2000',
                'description':
                'Bottom boundary layer model for the bed shear stress driving resuspension.',
                'level': CONFIG_LEVEL_BASIC
            }})

        # Which SG2000 shear stress drives resuspension: 'combined' (peak
        # wave-current, tau=rho*ustarcw^2, recommended for initiation of motion),
        # 'mean' (time-mean current, ustarc), or 'wave' (max wave, ustarwm).
        self._add_config({
            'vertical_mixing:bbl_stress': {
                'type': 'enum',
                'enum': ['combined', 'mean', 'wave'],
                'default': 'combined',
                'description':
                'Which SG2000 shear stress is compared against tau_crit for resuspension.',
                'level': CONFIG_LEVEL_BASIC
            }})

        # SG2000/log-law validity guard: the near-bed reference velocity must lie
        # within a resolvable bottom boundary layer. Where the deepest velocity cell
        # centre sits higher than this above the bed (coarse deep-ocean grid, e.g. a
        # ~50 m near-bottom layer), the log-law is invalid and the stress reverts to
        # the MITgcm-style quadratic drag (see _bottom_stress_sg2000).
        self._add_config({
            'vertical_mixing:bbl_max_ref_height': {
                'type': 'float',
                'default': 10.0,
                'min': 0.1,
                'max': 1000.0,
                'units': 'm',
                'description':
                'Max height of the near-bed reference velocity for which the BBL '
                'log-law / SG2000 is applied; above it the bed stress uses a '
                'quadratic drag (MITgcm bottom BC).',
                'level': CONFIG_LEVEL_BASIC
            }})
        self._add_config({
            'vertical_mixing:bottom_drag_coefficient': {
                'type': 'float',
                'default': 0.0021,
                'min': 0,
                'max': 0.1,
                'units': '1',
                'description':
                'Quadratic bottom drag coefficient c_d used where the BBL is not '
                'grid-resolved (and by the legacy scheme).',
                'level': CONFIG_LEVEL_BASIC
            }})

        # Bed/fluid properties that set the SG2000 bottom stress. These describe the
        # ambient seabed and water, NOT the drifting tracer: the bed shear stress is a
        # property of the flow + bed, so it must not depend on an individual particle's
        # grain size or density. (Per-particle grain_diameter/rho_s/tau_crit still drive
        # that particle's own settling, resuspension height, and resuspension threshold.)
        self._add_config({
            'vertical_mixing:bed_median_grain_size': {
                'type': 'float',
                'default': 31e-6,
                'min': 1e-7,
                'max': 1e-2,
                'units': 'm',
                'description':
                'Median grain size of the ambient seabed, used for SG2000 skin '
                'friction, ripple geometry and hydraulic roughness.',
                'level': CONFIG_LEVEL_BASIC
            }})
        self._add_config({
            'vertical_mixing:bed_sediment_density': {
                'type': 'float',
                'default': 2650.0,
                'min': 1000.0,
                'max': 5000.0,
                'units': 'kgm-3',
                'description':
                'Density of the ambient seabed sediment, used for the SG2000 density '
                'ratio s in skin friction / ripple roughness.',
                'level': CONFIG_LEVEL_BASIC
            }})

        # How the near-bed REFERENCE velocity (input to every BBL stress scheme) is
        # obtained for settled elements. The reader call to get it dominates runtime
        # at large settled counts, so this trades cost vs faithfulness:
        #   'profile'      - interpolate the full vertical profile, then pick the
        #                    deepest velocity cell above the seafloor (original; exact
        #                    but O(n_settled x n_levels) interpolation every step).
        #   'deepest_cell' - same deepest-cell value, but obtained with a SINGLE-level
        #                    reader call at the (cached, static) cell depth -> O(n_settled),
        #                    ~n_levels x cheaper, reproduces 'profile' to interpolation
        #                    round-off. Recommended for large-domain settling runs.
        #   'fixed_height' - near-bed velocity at a FIXED height above the bed
        #                    (bbl_fixed_ref_height); single-level call. Cheapest and
        #                    uses a consistent reference height, but changes the value
        #                    vs 'profile' (different reference height per cell).
        self._add_config({
            'vertical_mixing:bbl_ref_velocity': {
                'type': 'enum',
                'enum': ['profile', 'deepest_cell', 'fixed_height'],
                # Default 'deepest_cell': bit-identical to 'profile' (validated max|Δ|=0
                # Pa in bed stress) but a single-level reader call instead of a full
                # vertical-profile interpolation -- ~23% cheaper resuspension / bottom-
                # stress per step at large settled counts (see BOTTLENECK_LARGEDOMAIN.md
                # deepest_cell FOLLOW-UP). Set 'profile' to restore the exact original
                # interpolation path.
                'default': 'deepest_cell',
                'description':
                'How the near-bed reference velocity for the BBL stress is obtained.',
                'level': CONFIG_LEVEL_ADVANCED
            }})
        self._add_config({
            'vertical_mixing:bbl_fixed_ref_height': {
                'type': 'float',
                'default': 5.0,
                'min': 0.1,
                'max': 100.0,
                'units': 'm',
                'description':
                "Height above the bed at which the near-bed reference velocity is "
                "read when bbl_ref_velocity == 'fixed_height'.",
                'level': CONFIG_LEVEL_ADVANCED
            }})

        # ----- Dynamic critical shear stress (resuspension threshold) -----
        # How tau_crit is obtained per element (Update 1). Production default 'auto'
        # routes by sed_class; 'constant' keeps the legacy prescribed per-element
        # tau_crit (set this to reproduce old runs). The dynamic modes derive
        # tau_crit from particle properties + the ambient fluid:
        #   'shields'  - all elements: Soulsby-Whitehouse Shields from grain d, rho_s.
        #   'cohesive' - all elements: fractal floc strength (phi, Df, d_floc, rho_s).
        #   'auto'     - per element by sed_class (0->shields, 1->cohesive floc).
        #   'mixed'    - bed-composition blend Pc(bed_mud_fraction) of shields & floc.
        self._add_config({
            'vertical_mixing:tau_crit_mode': {
                'type': 'enum',
                'enum': ['constant', 'shields', 'cohesive', 'auto', 'mixed'],
                'default': 'auto',
                'description':
                'How the resuspension threshold tau_crit is computed per element. '
                "Production default 'auto' routes each element by sed_class "
                "(0->Shields, 1->cohesive floc strength); use 'constant' to "
                'reproduce legacy runs with a prescribed per-element tau_crit.',
                'level': CONFIG_LEVEL_BASIC
            }})
        self._add_config({
            'vertical_mixing:cohesive_size_cutoff': {
                'type': 'float',
                'default': 63e-6,
                'min': 0,
                'max': 1e-3,
                'units': 'm',
                'description':
                'Grain size below which an element left at sed_class=0 is treated '
                'as cohesive (sed_class=1) by the dynamic tau_crit closure.',
                'level': CONFIG_LEVEL_ADVANCED
            }})
        self._add_config({
            'vertical_mixing:cohesive_strength_coeff': {
                'type': 'float',
                'default': 100.0,
                'min': 0,
                'max': 1e6,
                'units': '1',
                'description':
                'Calibration coefficient c_str in the fractal cohesive strength '
                'tau = c_str*g*(rho_s-rho_f)*d_floc*phi^(2/(3-Df)). The buoyant '
                'stress scale g*(rho_s-rho_f)*d_floc is tiny for fine flocs, so '
                'c_str is O(100); the default gives ~0.04 Pa for the default fresh '
                'floc (d=4um, phi0=0.1, Df=2), comparable to the legacy tau_crit. '
                'Calibrate to your seeded floc properties / measured erosion threshold.',
                'level': CONFIG_LEVEL_ADVANCED
            }})
        self._add_config({
            'vertical_mixing:consolidation_timescale': {
                'type': 'float',
                'default': 86400.0,
                'min': 1.0,
                'max': 1e9,
                'units': 's',
                'description':
                'Timescale over which a settled cohesive deposit consolidates '
                '(phi grows phi0 -> consolidated_solids_fraction), raising tau_crit.',
                'level': CONFIG_LEVEL_ADVANCED
            }})
        self._add_config({
            'vertical_mixing:consolidated_solids_fraction': {
                'type': 'float',
                'default': 0.4,
                'min': 0,
                'max': 1,
                'units': '1',
                'description':
                'Packed solids volume fraction phi_max approached by a fully '
                'consolidated cohesive deposit.',
                'level': CONFIG_LEVEL_ADVANCED
            }})
        self._add_config({
            'vertical_mixing:mud_fraction_lower': {
                'type': 'float',
                'default': 0.05,
                'min': 0,
                'max': 1,
                'units': '1',
                'description':
                'Bed mud fraction at which cohesive behaviour begins (Pc=0 below) '
                'for the mixed tau_crit mode.',
                'level': CONFIG_LEVEL_ADVANCED
            }})
        self._add_config({
            'vertical_mixing:mud_fraction_upper': {
                'type': 'float',
                'default': 0.35,
                'min': 0,
                'max': 1,
                'units': '1',
                'description':
                'Bed mud fraction at which behaviour is fully cohesive (Pc=1 above) '
                'for the mixed tau_crit mode (Yao 2022 ~0.35 silt; Jacobs 2011 clay).',
                'level': CONFIG_LEVEL_ADVANCED
            }})
        self._add_config({
            'vertical_mixing:bed_mud_fraction': {
                'type': 'float',
                'default': 0.5,
                'min': 0,
                'max': 1,
                'units': '1',
                'description':
                'Mud (cohesive) fraction of the ambient seabed, used by the mixed '
                'tau_crit mode Pc blend.',
                'level': CONFIG_LEVEL_BASIC
            }})

        # ----- Resuspension height (where a resuspended element is placed) -----
        # How the height above the bed is set when an element resuspends (Update 2):
        #   'legacy'    - the original lift/gravity/drag ODE (vectorized solve_ivp),
        #                 integrated for a fixed 120 s (reproduces existing runs).
        #   'analytic'  - the closed-form force-balance EQUILIBRIUM of that same ODE
        #                 (no solver, ~zero cost). NB the ODE is strongly overdamped
        #                 and does not reach this equilibrium within the 120 s of
        #                 'legacy', so for fine sediment 'analytic' gives larger
        #                 heights than 'legacy' (same u* scaling, different magnitude).
        #                 CAVEAT: for fine/cohesive sediment the ODE equilibrium is
        #                 only reached over hours-to-days, so placing a particle there
        #                 in one minute-to-hour step over-resuspends; 'analytic' is
        #                 only timestep-consistent for fast-settling (coarse) grains.
        #                 Prefer 'turbulent' for fine/cohesive material.
        #   'reference' - place at a small reference height z_a and let the model's
        #                 vertical_mixing() carry it (turbulence-settling balance).
        #   'turbulent' - stochastic draw from the near-bed Rouse profile, with
        #                 near-bed turbulence parameterized from u* (hybrid: this
        #                 seeds the height, vertical_mixing() then evolves it).
        self._add_config({
            'vertical_mixing:resuspension_height_mode': {
                'type': 'enum',
                'enum': ['legacy', 'analytic', 'reference', 'turbulent'],
                'default': 'turbulent',
                'description':
                'How the above-bed height of a resuspended element is determined. '
                "Production default 'turbulent' (physical near-bed Rouse draw, "
                "timestep-consistent); use 'legacy' to reproduce old ODE runs.",
                'level': CONFIG_LEVEL_BASIC
            }})
        self._add_config({
            'vertical_mixing:resuspension_reference_height': {
                'type': 'float',
                'default': 0.01,
                'min': 1e-4,
                'max': 100.0,
                'units': 'm',
                'description':
                'Reference height z_a above the bed at which a resuspended element '
                'is placed (reference mode) or from which the Rouse draw starts '
                '(turbulent mode). Capped at 0.1*water_depth.',
                'level': CONFIG_LEVEL_BASIC
            }})
        self._add_config({
            'vertical_mixing:resuspension_seed_layer': {
                'type': 'float',
                'default': 5.0,
                'min': 0.0,
                'max': 1000.0,
                'units': 'm',
                'description':
                'Upper bound on the thickness above z_a over which the turbulent-mode '
                'Rouse draw can place a resuspended element in one step; the actual '
                'cap is min(this, kappa*u*dt), i.e. never more than turbulent '
                'diffusion can carry the particle in one time step. '
                'vertical_mixing() evolves it further afterwards.',
                'level': CONFIG_LEVEL_ADVANCED
            }})
        self._add_config({
            'vertical_mixing:rouse_beta': {
                'type': 'float',
                'default': 1.0,
                'min': 0.1,
                'max': 5.0,
                'units': '1',
                'description':
                'Sediment/momentum turbulent diffusivity ratio (inverse turbulent '
                'Schmidt number) beta in the Rouse number P=w_s/(beta*kappa*u*).',
                'level': CONFIG_LEVEL_ADVANCED
            }})
        self._add_config({
            'vertical_mixing:rouse_suspension_number': {
                'type': 'float',
                'default': 2.5,
                'min': 0.5,
                'max': 10.0,
                'units': '1',
                'description':
                'Rouse number P above which an element cannot be suspended and is '
                'kept in the near-bed (bedload) layer at z_a (turbulent mode).',
                'level': CONFIG_LEVEL_ADVANCED
            }})

        # By default, sediments do not strand towards coastline
        # TODO: A more sophisticated stranding algorithm is needed
        self._set_config_default('general:coastline_action', 'previous')

        # Vertical mixing is enabled as default
        self._set_config_default('drift:vertical_mixing', True)

        self.use_parallel = True

    def update(self):
        """Update positions and properties of sediment particles.
        """
        self.resuspension()
        # Advecting here all elements, but want to soon add
        # possibility of not moving settled elements, until
        # they are resuspended. May then need to send a boolean
        # array to advection methods below
        # -- UPDATE: Lana/Saulo, we have added this functionality somewhat,
        # the code still seems to waste time on the calculation maybe.
        self.advect_ocean_current()

        self.vertical_advection()

        self.advect_wind()  # Wind shear in upper 10cm of ocean

        self.stokes_drift()

        if self.get_config('drift:vertical_mixing') is False:
            self.vertical_buoyancy()
            # self.vertical_advection()
        else:
            self.vertical_mixing() # including buoyancy and settling


        self.deactivate_elements_outofbounds()       

        # upwards_moving_particles = self.elements.counter == 1
        # Restore downwards velocity to resuspended particles
        # self.elements.terminal_velocity[upwards_moving_particles] = self.elements.terminal_velocity_default[upwards_moving_particles]
        # self.elements.counter[upwards_moving_particles] = 0

        # self.resuspension()
        self.deactivate_elements(self.elements.beached == 1, reason='beached')
        self.remove_deactivated_elements()

    def update_terminal_velocity(self, Tprofiles=None, Sprofiles=None, z_index=None):
        """Per-element settling (terminal) velocity from the configured closure.

        Overrides the OceanDrift stub, which the base model already calls every
        step inside vertical_mixing() (and once before the run). The closure
        (config 'vertical_mixing:settling_model') maps each element's innate
        properties (grain_diameter, rho_s, fractal_dim, d0, corey_shape_factor)
        and the *local* fluid (rho_f from T/S, kinematic viscosity from the
        molecular viscosity config) to a settling speed via models/settling.py;
        terminal_velocity is stored negative (downward).

        Static vs dynamic (config 'vertical_mixing:settling_dynamic'):
          * False (default) - compute only for elements not yet set (vel_set==0),
            i.e. once at first activation, then hold fixed. Cheap, seed-time static.
          * True - recompute for all active elements every step (rho_f, nu follow
            the local T/S), giving a dynamic terminal velocity with no other change.

        The legacy 'prescribed' model is a no-op: terminal_velocity keeps whatever
        was seeded / set by SedimentElement.move_elements (back-compat).
        """
        model = self.get_config('vertical_mixing:settling_model')
        if model == 'prescribed':
            return

        if self.get_config('vertical_mixing:settling_dynamic'):
            idx = np.ones(self.elements.z.shape, dtype=bool)
        else:
            idx = self.elements.vel_set == 0
        if not np.any(idx):
            return

        # Local fluid: density from in-situ T/S, kinematic viscosity from the
        # configured molecular (dynamic) viscosity. Same provenance as the SG2000
        # bed stress and the dynamic tau_crit closure.
        temp = self.environment['sea_water_temperature'][idx]
        sal = self.environment['sea_water_salinity'][idx]
        rho_f = np.broadcast_to(
            np.asarray(self.sea_water_density(temp, sal), dtype=float),
            (int(np.count_nonzero(idx)),)).astype(float)
        nu = self.get_config('environment:molecular_viscosity') / rho_f  # m^2/s

        w = settling.settling_velocity(
            model,
            d=self.elements.grain_diameter[idx].astype(float),
            rho_s=self.elements.rho_s[idx].astype(float),
            rho_f=rho_f, nu=nu,
            corey_shape_factor=self.elements.corey_shape_factor[idx].astype(float),
            d0=self.elements.d0[idx].astype(float),
            fractal_dim=self.elements.fractal_dim[idx].astype(float))

        self.elements.terminal_velocity[idx] = -w   # negative == sinking
        self.elements.vel_set[idx] = 1

    def deactivate_elements_outofbounds(self):
        # This only works if the first reader you passed to the model
        # has the correct geographic bounds
        reader_names = list(self.env.readers.keys())
        reader_name = reader_names[0]
        lon_min = self.env.readers[reader_name].xmin
        lon_max = self.env.readers[reader_name].xmax
        lat_min = self.env.readers[reader_name].ymin
        lat_max = self.env.readers[reader_name].ymax
        
        lons, lats = self.elements.lon, self.elements.lat
        out_of_bounds = (lons < lon_min) | (lons > lon_max) | (lats < lat_min) | (lats > lat_max)
        self.deactivate_elements(out_of_bounds, reason='out of bounds')
        # self.remove_deactivated_elements()

    def bottom_interaction(self, seafloor_depth):
        """Sub method of vertical_mixing and vertical_buoyancy, determines settling"""
        # Elements at or below seafloor are settled, by setting
        # self.elements.moving to 0.
        # These elements will not move until eventual later resuspension.
        settling = np.logical_and(self.elements.z <= seafloor_depth, self.elements.moving==1)
        if np.sum(settling) > 0:
            logger.debug('Settling %s elements at seafloor' % np.sum(settling))
            self.elements.moving[settling] = 0
            self.elements.settled[settling] = 1
            # Freshly deposited: restart the consolidation clock so the cohesive
            # tau_crit begins from its as-deposited (phi0) value.
            self.elements.time_since_settled[settling] = 0.0
            #  self.get_distance_cell_center_above(settling)  # no longer needed here

    def _critical_shear_stress(self, idxs):
        """Resuspension threshold tau_crit [Pa] for the masked elements (Update 1).

        Returned as a full-size array. The threshold is a function of each
        element's *innate* properties (grain size/density, sediment class, floc
        fractal state) and the *ambient fluid* (density from local T/S, molecular
        viscosity) -- never of the bed shear stress, which is computed separately.

        Modes (config 'vertical_mixing:tau_crit_mode'):
          'constant' - legacy prescribed per-element tau_crit (back-compat default).
          'shields'  - Soulsby-Whitehouse Shields for all elements.
          'cohesive' - fractal floc strength for all elements.
          'auto'     - per element by sed_class (0 -> shields, 1 -> floc); elements
                       left at sed_class=0 with grain size below
                       'cohesive_size_cutoff' are treated as cohesive.
          'mixed'    - bed-composition blend Pc(bed_mud_fraction) of shields & floc.
        """
        tau = self.elements.tau_crit.astype(float).copy()  # constant / fallback
        mode = self.get_config('vertical_mixing:tau_crit_mode')
        if mode == 'constant':
            return tau

        fidx = np.where(idxs)[0]
        if fidx.size == 0:
            return tau

        # Ambient fluid on the subset (same provenance as the SG2000 bed stress).
        temp = self.environment['sea_water_temperature'][fidx]
        sal = self.environment['sea_water_salinity'][fidx]
        rho_f = np.broadcast_to(
            np.asarray(self.sea_water_density(temp, sal), dtype=float),
            (fidx.size,)).astype(float)
        nu_dyn = self.get_config('environment:molecular_viscosity')  # kg/m/s

        # Element (innate) properties on the subset.
        d = self.elements.grain_diameter[fidx].astype(float)
        rho_s = self.elements.rho_s[fidx].astype(float)

        # Non-cohesive (Shields) threshold.
        tau_nc = floc_strength.tau_crit_shields(d, rho_s, rho_f, nu_dyn)

        # Cohesive (fractal floc) threshold, with consolidation growth of phi.
        Df = self.elements.fractal_dim[fidx].astype(float)
        phi0 = self.elements.phi0[fidx].astype(float)
        t_settled = self.elements.time_since_settled[fidx].astype(float)
        phi = floc_strength.consolidate_phi(
            phi0, t_settled,
            self.get_config('vertical_mixing:consolidation_timescale'),
            self.get_config('vertical_mixing:consolidated_solids_fraction'))
        tau_coh = floc_strength.tau_crit_floc(
            phi, Df, d, rho_s, rho_f,
            c_str=self.get_config('vertical_mixing:cohesive_strength_coeff'))

        if mode == 'shields':
            tau[fidx] = tau_nc
        elif mode == 'cohesive':
            tau[fidx] = tau_coh
        elif mode == 'auto':
            cutoff = self.get_config('vertical_mixing:cohesive_size_cutoff')
            is_cohesive = (self.elements.sed_class[fidx] == 1) | (d < cutoff)
            tau[fidx] = np.where(is_cohesive, tau_coh, tau_nc)
        elif mode == 'mixed':
            fmud = self.get_config('vertical_mixing:bed_mud_fraction')
            lo = self.get_config('vertical_mixing:mud_fraction_lower')
            hi = self.get_config('vertical_mixing:mud_fraction_upper')
            Pc = np.clip((fmud - lo) / max(hi - lo, 1e-6), 0.0, 1.0)
            tau[fidx] = (1.0 - Pc) * tau_nc + Pc * tau_coh

        return tau

    def resuspension(self):
        """Resuspending elements when critical shear stress is exceed"""
        ## Old resuspension condition
        # threshold = self.get_config('vertical_mixing:resuspension_threshold')
        # resuspending = np.logical_and(self.current_speed() > threshold, self.elements.moving==0)

        settled = np.logical_and(self.elements.moving == 0, self.elements.settled == 1)
        if np.sum(settled) > 0:
            # Advance the consolidation clock for elements resting on the bed.
            # (Used by the cohesive branch of the dynamic tau_crit closure; a no-op
            # for tau_crit_mode='constant'/'shields'.)
            self.elements.time_since_settled[settled] = (
                self.elements.time_since_settled[settled] + self.time_step.total_seconds())
            # Dynamic threshold: a function of each element's properties and the
            # ambient fluid (Update 1). 'constant' mode returns the prescribed
            # per-element tau_crit, i.e. the legacy behaviour.
            threshold = self._critical_shear_stress(settled)
            bottom_stress = self.calc_bottom_stress(settled)
            resuspending = np.logical_and(bottom_stress > threshold, settled)
        else:
            return
            #  resuspending = settled
        
        if np.sum(resuspending) > 0:
            # Keep track of how many times particle has been resuspended
            self.elements.times_resuspended[resuspending] = self.elements.times_resuspended[resuspending] + 1
            # Allow moving again
            self.elements.moving[resuspending] = 1
            # Since not at bottom anymore, set distance to nearest cell center above to 99999
            # self.elements.dz_bot[resuspending] = 99999.0  # no longer needed here
            # Give particle upwards velocity
            # self.elements.terminal_velocity[resuspending] = self.calc_upward_resuspension_velocity(bottom_stress, resuspending)
            # Keep track of number of time steps particle has been in resuspension for (however, it is always 1 dt)
            # self.elements.counter[resuspending] = 1  # should be improved - a nominal dt is too long for being resuspended, maybe as a var called in_resuspension that sends these elements to a small loop like v-mix
            # Calculate the particle z position from height above bottom following
            # the configured resuspension-height scheme (Update 2).
            height_above_bottom = self.calc_resuspension_height(bottom_stress, resuspending)
            # print(np.max(height_above_bottom))
            new_z = self.elements.z[resuspending] + height_above_bottom
            new_z[new_z > 0] = 0
            # Update the z position of the particle
            self.elements.z[resuspending] = new_z
            # Switch particles that are resuspending out of settled mode
            self.elements.settled[resuspending] = 0
            # Reset the consolidation clock: a resuspended element re-deposits fresh.
            self.elements.time_since_settled[resuspending] = 0.0
            # keep track of the latest resuspension height:
            self.elements.latest_resuspension_height[resuspending] = height_above_bottom

    def _near_bottom_environment(self, idxs):
        """Vectorized extraction of the near-bed environment for the masked
        elements `idxs`. This factors the MITgcm-grid / input logic shared by all
        BBL schemes: a vectorized `searchsorted` locates the deepest velocity cell
        above the local seafloor, from which the near-bed current, vertical
        diffusivity and layer thickness `dz_bot` are read. Surface wave fields and
        water depth are taken from `self.environment` (already computed this step,
        so no extra reader call). Elements settled below the shallowest grid cell
        are flagged `beached`.

        Returns a dict of subset-length (n_active) arrays, or None if empty.
        """
        floor_idx = np.where(idxs)[0]
        if floor_idx.size == 0:
            return None

        lons = self.elements.lon[idxs]
        lats = self.elements.lat[idxs]
        atimes = self.time
        n_active = floor_idx.size
        dz_bot_min = 0.1
        # Seafloor depth (negative z) is already in self.environment from this step.
        el_seafloor_depth_z = -1 * self.environment['sea_floor_depth_below_sea_level'][idxs]

        variables = ['x_sea_water_velocity', 'y_sea_water_velocity',
                     'ocean_vertical_diffusivity', 'sea_floor_depth_below_sea_level']
        scheme = self.get_config('vertical_mixing:bbl_ref_velocity')

        # ---- 'fixed_height' (A): near-bed velocity at a fixed height above the bed ----
        if scheme == 'fixed_height':
            h_ref = float(self.get_config('vertical_mixing:bbl_fixed_ref_height'))
            # reference z: h_ref above the bed, but kept below the surface
            z_ref = np.minimum(el_seafloor_depth_z + h_ref, -dz_bot_min)
            env_var, _, _ = self.env.get_environment(
                variables, atimes, lons, lats, z_ref, None)   # profiles=None -> single level
            return dict(floor_idx=floor_idx, n_active=n_active,
                        u=np.asarray(env_var['x_sea_water_velocity'], dtype=float),
                        v=np.asarray(env_var['y_sea_water_velocity'], dtype=float),
                        visc=np.asarray(env_var['ocean_vertical_diffusivity'], dtype=float),
                        dz_bot=np.full(n_active, h_ref),
                        depth=-el_seafloor_depth_z)

        # 'profile' and 'deepest_cell' both locate the deepest velocity cell above the
        # seafloor on the reader's (static) vertical grid `zgrid`.
        profiles = ['x_sea_water_velocity', 'y_sea_water_velocity', 'ocean_vertical_diffusivity']
        need_profile_call = (scheme == 'profile') or (getattr(self, '_bbl_zgrid', None) is None)
        if need_profile_call:
            # Full-profile reader call: the original path, and also how we capture and
            # cache the static reader z-grid on the first 'deepest_cell' step.
            env_var, env_profs, amiss = self.env.get_environment(
                variables, atimes, self.elements.lon[idxs], lats, self.elements.z[idxs], profiles)
            self._bbl_zgrid = np.asarray(env_profs['z'])

        zgrid = self._bbl_zgrid
        search_indices = np.searchsorted(-zgrid, -el_seafloor_depth_z, side='right')
        selected_idx = np.where(search_indices > 0, search_indices - 1, 0)
        dz_bot = np.maximum(zgrid[selected_idx] - el_seafloor_depth_z, dz_bot_min)
        fallback = (search_indices == 0)
        if np.any(fallback):
            self.elements.beached[floor_idx[fallback]] = 1

        if scheme == 'profile':
            particle_range = np.arange(n_active)
            u = env_profs['x_sea_water_velocity'][selected_idx, particle_range]
            v = env_profs['y_sea_water_velocity'][selected_idx, particle_range]
            visc = env_profs['ocean_vertical_diffusivity'][selected_idx, particle_range]
        else:
            # ---- 'deepest_cell' (A'): single-level reader call at the cell depth ----
            # zgrid[selected_idx] is the deepest velocity-cell depth above each
            # element's seafloor; reading there (profiles=None) reproduces the
            # 'profile' pick without interpolating the whole column.
            z_ref = zgrid[selected_idx]
            env_var, _, _ = self.env.get_environment(
                variables, atimes, lons, lats, z_ref, None)
            u = np.asarray(env_var['x_sea_water_velocity'], dtype=float)
            v = np.asarray(env_var['y_sea_water_velocity'], dtype=float)
            visc = np.asarray(env_var['ocean_vertical_diffusivity'], dtype=float)

        return dict(floor_idx=floor_idx, n_active=n_active,
                    u=u, v=v, visc=visc, dz_bot=dz_bot,
                    depth=-el_seafloor_depth_z)  # positive-down water depth [m]

    @staticmethod
    def _wave_number(omega, h, g=9.81, n_iter=8):
        """Linear-wave dispersion: solve omega^2 = g k tanh(k h) for k.
        Vectorized Newton iteration from the deep-water guess."""
        omega = np.asarray(omega, dtype=float)
        h = np.asarray(h, dtype=float)
        k = np.maximum(omega ** 2 / g, 1e-6)
        for _ in range(n_iter):
            th = np.tanh(np.minimum(k * h, 50.0))
            f = g * k * th - omega ** 2
            df = g * th + g * k * h * (1 - th ** 2)
            k = np.maximum(k - f / df, 1e-8)
        return k

    def _bottom_stress_sg2000(self, env):
        """Styles & Glenn (2000) bed shear stress [Pa] for the near-bed `env`.

        The bed stress is computed by one of three routes per element, chosen by
        whether the near-bed reference velocity height (`dz_bot`, distance from the
        seafloor to the centre of the deepest velocity cell) lies within a
        resolvable bottom boundary layer:

          1. dz_bot > `vertical_mixing:bbl_max_ref_height`  -> the reference current
             sits ABOVE the log/BBL layer (coarse deep-ocean grid, e.g. a ~50 m thick
             near-bottom cell). The log-law / SG2000 inversion is not valid there, so
             the stress reverts to the MITgcm-style quadratic drag tau = rho*c_d*|u|^2
             (same rationale and c_d as the legacy scheme / MITgcm bottom BC).
          2. dz_bot <= cap and no waves at the bed -> current-only law-of-the-wall.
          3. dz_bot <= cap and waves reach the bed (Ub>1mm/s) -> full SG2000 solve.

        So SG2000/log-law is applied only where the grid resolves the BBL (typically
        the shelf), while deep coarse-grid cells use the drag law. Fully vectorized
        (no per-element Python loops); SG2000 internals are cgs.
        """
        g = 9.81
        kappa = 0.4
        fidx = env['floor_idx']
        n = env['n_active']
        bbl_max_ref_height = self.get_config('vertical_mixing:bbl_max_ref_height')
        c_d = self.get_config('vertical_mixing:bottom_drag_coefficient')

        Ur = np.maximum(np.sqrt(env['u'] ** 2 + env['v'] ** 2), 1e-6)  # m/s
        zr = np.maximum(env['dz_bot'], 1e-3)                           # m
        depth = np.maximum(env['depth'], 1e-3)                         # m, positive down

        # Bed/fluid properties (environment + config only) — the bed shear stress must
        # NOT depend on the individual drifting particle. rho_f from the local T/S;
        # bed median grain size, bed sediment density and molecular viscosity from config.
        temp = self.environment['sea_water_temperature'][fidx]
        sal = self.environment['sea_water_salinity'][fidx]
        rho_f = np.broadcast_to(
            np.asarray(self.sea_water_density(temp, sal), dtype=float), (n,)).astype(float)
        d_med = np.full(n, self.get_config('vertical_mixing:bed_median_grain_size'))  # m
        rho_bed = self.get_config('vertical_mixing:bed_sediment_density')             # kg/m3
        nu_dyn = self.get_config('environment:molecular_viscosity')                   # kg/m/s

        # wave forcing (surface fields, from self.environment)
        Hs = np.asarray(self.environment['sea_surface_wave_significant_height'])[fidx].astype(float)
        Tp = np.asarray(self.environment[
            'sea_surface_wave_period_at_variance_spectral_density_maximum'])[fidx].astype(float)

        # near-bed wave orbital velocity / excursion from linear-wave theory.
        # In deep water (e.g. the San Pedro basin) waves do not reach the bed and Ub->0;
        # those cells fall back to the current-only law-of-the-wall. When there is no
        # wave forcing at all (Hs==0, the typical MITgcm-only run) the whole dispersion
        # computation is skipped so the cost matches the current-only legacy estimate.
        a_orb = np.zeros(n)
        Ub_w = np.zeros(n)
        if np.any((Hs > 1e-4) & (Tp > 1e-2)):
            Tp_ok = Tp > 1e-2
            omega = np.where(Tp_ok, 2.0 * np.pi / np.where(Tp_ok, Tp, 1.0), 0.0)
            kwave = self._wave_number(omega, depth, g)
            sinhkh = np.maximum(np.sinh(np.minimum(kwave * depth, 50.0)), 1e-6)
            a_orb = np.where(Tp_ok & (Hs > 0), (Hs / 2.0) / sinhkh, 0.0)  # excursion amp [m]
            Ub_w = a_orb * omega                                          # orbital vel [m/s]
        wave_active = (Ub_w > 1e-3)   # >~1 mm/s near-bed orbital velocity

        # Is the reference velocity within a resolvable bottom boundary layer?
        bbl_valid = zr <= bbl_max_ref_height

        ustarcw = np.zeros(n)
        ustarc = np.zeros(n)
        ustarwm = np.zeros(n)

        # ---- (1) reference height above the BBL -> MITgcm quadratic drag ----
        # tau = rho*c_d*|u|^2  i.e. u* = sqrt(c_d)*|u| (legacy / MITgcm bottom BC).
        drag = ~bbl_valid
        if np.any(drag):
            us = np.sqrt(c_d) * Ur[drag]
            ustarc[drag] = us
            ustarcw[drag] = us

        # ---- (2) current-only law-of-the-wall where the BBL is resolved ----
        loglaw = bbl_valid & (~wave_active)
        if np.any(loglaw):
            kbr_def_m = 0.03  # 3 cm default bed roughness (SG2000 no-motion limit)
            z0 = (d_med[loglaw] + kbr_def_m) / 30.0
            ratio = np.maximum(zr[loglaw] / np.maximum(z0, 1e-9), 1.0001)
            us = kappa * Ur[loglaw] / np.log(ratio)
            ustarc[loglaw] = us
            ustarcw[loglaw] = us

        # ---- (3) full SG2000 where the BBL is resolved and waves reach the bed ----
        sg = bbl_valid & wave_active
        if np.any(sg):
            wa = sg
            # wave-current angle from stokes-drift vs near-bed current direction
            stx = np.asarray(self.environment[
                'sea_surface_wave_stokes_drift_x_velocity'])[fidx][wa].astype(float)
            sty = np.asarray(self.environment[
                'sea_surface_wave_stokes_drift_y_velocity'])[fidx][wa].astype(float)
            cur_dir = np.arctan2(env['v'][wa], env['u'][wa])
            wav_dir = np.arctan2(sty, stx)
            dtheta = np.arctan2(np.sin(wav_dir - cur_dir), np.cos(wav_dir - cur_dir))
            # undefined direction (no stokes info) -> assume aligned
            dtheta = np.where(np.hypot(stx, sty) < 1e-9, 0.0, dtheta)
            deg = np.degrees(np.abs(dtheta))

            # bed/fluid properties on the wave-active subset (config-based, per-cell rho_f)
            s_wa = rho_bed / rho_f[wa]
            nu_wa = nu_dyn / rho_f[wa]

            out = sg2000_solve(
                Ub=Ub_w[wa] * 100.0, Ab=a_orb[wa] * 100.0, Ur=Ur[wa] * 100.0,
                zr=zr[wa] * 100.0, deg=deg,
                d_median=d_med[wa] * 100.0, s=s_wa,
                nu=nu_wa * 1e4, g=g * 100.0)
            ustarcw[wa] = out['ustarcw'] / 100.0   # cm/s -> m/s
            ustarc[wa] = out['ustarc'] / 100.0
            ustarwm[wa] = out['ustarwm'] / 100.0

        which = self.get_config('vertical_mixing:bbl_stress')
        ustar = {'combined': ustarcw, 'mean': ustarc, 'wave': ustarwm}[which]
        return rho_f * ustar ** 2   # Pa

    def calc_bottom_stress(self, idxs, turnoff=False):
        """Bed shear stress [Pa] for the masked elements, returned in a full-size
        array (zero outside `idxs`).

        Scheme selected by config 'vertical_mixing:bbl_scheme':
          'sg2000' - Styles & Glenn (2000) combined wave-current BBL (current-only
                     law-of-the-wall where there are no waves). See
                     `_bottom_stress_sg2000`.
          'legacy' - the previous c_d drag + diffusivity estimate.
        Both share the vectorized near-bed/grid extraction in
        `_near_bottom_environment`.
        """
        if turnoff:
            return np.zeros(self.elements.z.shape)

        env = self._near_bottom_environment(idxs)
        if env is None:
            return np.zeros(self.elements.z.shape)

        scheme = self.get_config('vertical_mixing:bbl_scheme')
        if scheme == 'legacy':
            r_b = 0
            c_d = 0.0021
            crit_last_layer_thickness = 1.
            vel_mag = np.sqrt(env['u'] ** 2 + env['v'] ** 2)
            turb_contrib = np.zeros(env['n_active'])
            turb_layers = env['dz_bot'] >= crit_last_layer_thickness
            turb_contrib[turb_layers] = (2 * env['visc'][turb_layers]
                                         / env['dz_bot'][turb_layers] * vel_mag[turb_layers])
            mean_flow_contrib = (r_b + c_d * vel_mag) * vel_mag
            computed = self.sea_water_density() * (mean_flow_contrib + turb_contrib)
        else:
            computed = self._bottom_stress_sg2000(env)

        full_bottom_stress = np.zeros(self.elements.z.shape)
        full_bottom_stress[env['floor_idx']] = computed
        return full_bottom_stress

    # def calc_bottom_stress(self, idxs):
    #     """
    #     Compute bottom stress for a set of particles (or grid cells) identified by the boolean mask `idxs`.
    #     The returned array has the same shape as the full elements array. For indices not included in `idxs`,
    #     the bottom stress is set to NaN.
        
    #     This function uses a vectorized search (via np.searchsorted) to locate the appropriate vertical
    #     cell for each active particle, then computes the bottom stress with the same logic as before.
    #     """
    #     # 1. Extract data only for active particles (where idxs is True)
    #     floor_idx = np.where(idxs)[0]  # indices in the full array that are active
    #     if floor_idx.size == 0:
    #         # No active particles: return an array of zeros with the full shape.
    #         return np.zeros(self.elements.z.shape)
            
    #     zs   = self.elements.z[idxs]
    #     lons = self.elements.lon[idxs]
    #     lats = self.elements.lat[idxs]
    #     atimes = self.time

    #     variables = ['x_sea_water_velocity', 'y_sea_water_velocity',
    #                 'ocean_vertical_diffusivity', 'sea_floor_depth_below_sea_level']
    #     profiles = ['x_sea_water_velocity', 'y_sea_water_velocity', 'ocean_vertical_diffusivity']

    #     # Get environment data (both cell-centered and profile)
    #     env_var, env_profs, amiss = self.env.get_environment(variables, atimes, lons, lats, zs, profiles)

    #     # seafloor depth (convert sign so that depth is positive downward)
    #     # (Note: self.environment['sea_floor_depth_below_sea_level'] is assumed to be defined over the entire domain)
    #     el_seafloor_depth_z = -1 * self.environment['sea_floor_depth_below_sea_level'][idxs]

    #     # Get the vertical grid from the profile data.
    #     # Here we assume that zgrid is a 1D array (e.g., shape (n_levels,)) that is monotonic.
    #     zgrid = env_profs['z']  # e.g., descending (surface to seafloor)

    #     # 2. Initialize arrays for active particles
    #     n_active = floor_idx.size
    #     # For each active particle, we want to determine the local near-bottom cell.
    #     # Set minimum bottom layer thickness and a critical layer thickness.
    #     dz_bot_min = 0.1                # minimum allowed thickness (m)
    #     crit_last_layer_thickness = 1.  # threshold for turbulent contribution

    #     # Allocate arrays of length n_active for near-bottom quantities.
    #     # (They will be filled with vectorized operations below.)
    #     u_bot    = np.empty(n_active)
    #     v_bot    = np.empty(n_active)
    #     visc_bot = np.empty(n_active)
    #     dz_bot   = np.full(n_active, dz_bot_min)
    #     turb_contrib = np.zeros(n_active)

    #     # 3. Find the near-bottom level for each active particle
    #     #
    #     # For each active particle with seafloor depth `adep` (from el_seafloor_depth_z),
    #     # we need the index of the last grid cell in zgrid that is above the seafloor.
    #     #
    #     # The original code did:
    #     #     last_indz = np.argwhere(zgrid > adep)
    #     # and used last_indz[-1] if any exist; otherwise, fallback to index 0.
    #     #
    #     # We can do this vectorially if we assume that zgrid is monotonic. However,
    #     # since np.searchsorted requires an ascending array, we invert the sign.
    #     #
    #     # For each active particle, we compute:
    #     #    search_idx = np.searchsorted(-zgrid, -adep, side='right')
    #     # If search_idx > 0 then the desired index is search_idx - 1; otherwise, we use 0.
    #     search_indices = np.searchsorted(-zgrid, -el_seafloor_depth_z, side='right')
    #     selected_idx = np.where(search_indices > 0, search_indices - 1, 0)

    #     # 4. Extract near-bottom values from the profile arrays.
    #     # We assume that env_profs['x_sea_water_velocity'], etc., have shape (n_levels, n_active)
    #     particle_range = np.arange(n_active)
    #     u_bot    = env_profs['x_sea_water_velocity'][selected_idx, particle_range]
    #     v_bot    = env_profs['y_sea_water_velocity'][selected_idx, particle_range]
    #     visc_bot = env_profs['ocean_vertical_diffusivity'][selected_idx, particle_range]

    #     # Compute the local layer thickness for each active particle.
    #     dz_bot = zgrid[selected_idx] - el_seafloor_depth_z

    #     # Log a debug message for any particle that got the fallback index.
    #     fallback = (search_indices == 0)
    #     if np.any(fallback):
    #         for i in np.where(fallback)[0]:
    #             self.elements.beached[floor_idx[i]] = 1
    #             logger.debug('Element %s settled at seafloor shallower than shallowest grid cell' %
    #                         self.elements.ID[floor_idx[i]])

    #     # 5. Compute the bottom stress using the same formulas as before
    #     r_b = 0
    #     c_d = 0.0021  # drag coefficient (to be updated later for a log-wall model)
    #     _2KE = u_bot**2 + v_bot**2  # twice the kinetic energy (ignoring vertical velocity)
    #     vel_mag = np.sqrt(_2KE)

    #     # Turbulent contribution is computed only if the layer thickness is above a threshold.
    #     turb_layers = dz_bot >= crit_last_layer_thickness
    #     turb_contrib[turb_layers] = 2 * visc_bot[turb_layers] / dz_bot[turb_layers] * vel_mag[turb_layers]

    #     # Mean-flow contribution (as in the original code)
    #     mean_flow_contrib = (r_b + c_d * np.sqrt(_2KE)) * vel_mag
    #     bottom_stress_pseudo_energy_magnitude = mean_flow_contrib + turb_contrib

    #     # Multiply by the sea water density.
    #     computed_bottom_stress = self.sea_water_density() * bottom_stress_pseudo_energy_magnitude

    #     # 6. Place the computed values back into an array of full size
    #     full_bottom_stress = np.zeros(self.elements.z.shape)
    #     # For the indices in the full array corresponding to active particles, insert the computed values.
    #     full_bottom_stress[floor_idx] = computed_bottom_stress

    #     return full_bottom_stress

    def find_nearest(self, array, value):
        """
        Find the nearest value in an array to a given value.

        Parameters
        ----------
        array : numpy.ndarray
            The array in which to search for the nearest value.
        value : float
            The value to find the nearest to in the array.

        Returns
        -------
        tuple
            A tuple containing the nearest value and its index in the array.
        """
        idx = (np.abs(array - value)).argmin()
        return array[idx], idx


    def calc_resuspension_height(self, bottom_stress, resuspending):
        """Above-bed height [m] for each resuspending element (Update 2).

        Dispatches on config 'vertical_mixing:resuspension_height_mode'. Returns a
        subset-length array (one value per True entry in `resuspending`), matching
        the legacy `calc_resuspension_depth` contract used by `resuspension()`.

          'legacy'    - the original lift/gravity/drag ODE (vectorized solve_ivp),
                        integrated 120 s (reproduces existing runs).
          'analytic'  - the closed-form force-balance equilibrium of that same ODE
                        z_eq = ln2*u* * sqrt(3 C_l rho_f / (8 r rho_s g (1-rho_f/rho_s))).
                        The ODE is overdamped and does not reach z_eq within the 120 s
                        of 'legacy', so this gives larger heights for fine sediment
                        (same u* scaling, different magnitude) -- it is the equilibrium
                        the ODE asymptotes to, not a reproduction of 'legacy'.
          'reference' - a small reference height z_a (then vertical_mixing carries it).
          'turbulent' - stochastic near-bed Rouse draw with u*-parameterized
                        turbulence (then vertical_mixing carries it).
        """
        mode = self.get_config('vertical_mixing:resuspension_height_mode')
        if mode == 'legacy':
            return self.calc_resuspension_depth(bottom_stress, resuspending)

        # Common near-bed quantities on the resuspending subset.
        temp = self.environment['sea_water_temperature'][resuspending]
        sal = self.environment['sea_water_salinity'][resuspending]
        rho_f = np.asarray(self.sea_water_density(temp, sal), dtype=float)
        bot = np.maximum(bottom_stress[resuspending].astype(float), 0.0)
        u_star = np.sqrt(bot / rho_f)
        n = int(np.count_nonzero(resuspending))

        if mode == 'analytic':
            g = 9.81
            C_l = self.elements.C_l[resuspending].astype(float)
            rho_s = self.elements.rho_s[resuspending].astype(float)
            r = self.elements.grain_diameter[resuspending].astype(float) / 2.0
            buoy = np.maximum(1.0 - rho_f / rho_s, 1e-6)
            return np.log(2.0) * u_star * np.sqrt(
                3.0 * C_l * rho_f / (8.0 * r * rho_s * g * buoy))

        # reference / turbulent: need the local water depth and reference height.
        h = np.maximum(
            self.environment['sea_floor_depth_below_sea_level'][resuspending].astype(float),
            0.1)
        za = np.minimum(self.get_config('vertical_mixing:resuspension_reference_height'),
                        0.1 * h)

        if mode == 'reference':
            return za

        # turbulent: stochastic Rouse draw (global seeded RNG via np.random).
        # Cap the single-step seed at the distance turbulent diffusion can actually
        # carry a particle in one model time step, kappa*u*dt, so the scheme stays
        # consistent with the time step (never seeds higher than physically
        # reachable per step); vertical_mixing() then evolves it further. This
        # makes 'turbulent' timestep-consistent for any dt (see RESUSPENSION_MECHANICS_PLAN.md).
        w_s = np.abs(self.elements.terminal_velocity[resuspending].astype(float))
        dt = self.time_step.total_seconds()
        seed_layer = np.minimum(
            self.get_config('vertical_mixing:resuspension_seed_layer'),
            rouse.KAPPA * u_star * dt)
        U = np.random.uniform(size=n)
        return rouse.resuspension_height(
            u_star, w_s, h, za, seed_layer=seed_layer,
            beta=self.get_config('vertical_mixing:rouse_beta'),
            P_crit=self.get_config('vertical_mixing:rouse_suspension_number'),
            U=U)

    def calc_resuspension_depth(self, bottom_stress, resuspending, rtol=1e-1, atol=1e-5, dameth='BDF'):
        """Legacy resuspension-height ODE (the `legacy` height mode).

        Integrates the lift/gravity/Stokes-drag particle-motion ODE for 120 s and
        returns the above-bed height. Same physics as before, but the previous
        joblib-process / batched implementation has been replaced by a single
        vectorized `solve_ivp` call with an analytic SPARSE Jacobian:

          * one solver call (no joblib process spawn/pickle overhead, no batching);
          * the per-particle dynamics are decoupled, so the 2n x 2n Jacobian is
            sparse (3n non-zeros). Supplying it analytically avoids the O(n)
            finite-difference RHS evaluations per Jacobian and lets BDF use sparse
            linear algebra, so the cost scales ~O(n) instead of the dense O((2n)^3).

        Results match the old batched version to within the solver tolerance.
        """
        gravity = 9.81
        temp = self.environment['sea_water_temperature'][resuspending]
        sal = self.environment['sea_water_salinity'][resuspending]
        rho_f = np.asarray(self.sea_water_density(temp, sal), dtype=float)

        C_l = self.elements.C_l[resuspending].astype(float)
        rho_s = self.elements.rho_s[resuspending].astype(float)
        r = self.elements.grain_diameter[resuspending].astype(float) / 2.0
        mu = self.elements.viscosity_molecular[resuspending].astype(float)
        u_star = np.sqrt(np.maximum(bottom_stress[resuspending].astype(float), 0.0) / rho_f)

        n = int(np.count_nonzero(resuspending))
        if n == 0:
            return np.array([])

        params = (C_l, rho_f, rho_s, u_star, r, mu, np.full(n, gravity))

        # Per-particle Jacobian constants: lift = B/z^2, drag term = drag_c * w.
        B = (0.5 * 3.0 * C_l * rho_f / (4.0 * r * rho_s)) * (u_star ** 2) * (np.log(2.0) ** 2)
        drag_c = 9.0 * mu / (2.0 * r ** 2 * rho_s)
        jrows, jcols = _resuspension_ode_jacobian_pattern(n)
        ones_n = np.ones(n)

        def jac(t, y):
            z = y[:n]
            z_fixed = np.where(z == 0, 1e-8, z)
            d_dwdt_dz = -2.0 * B / z_fixed ** 3
            data = np.concatenate([ones_n, d_dwdt_dz, -drag_c])
            return sparse.csr_matrix((data, (jrows, jcols)), shape=(2 * n, 2 * n))

        z0 = np.full(n, 1e-3 + 1e-8)
        y0 = np.concatenate([z0, np.zeros(n)])   # state ordered [z block, w block]

        sol = solve_ivp(
            fun=lambda t, y: particle_motion_vectorized(t, y, params),
            t_span=(0, 120), y0=y0, method=dameth, jac=jac,
            rtol=rtol, atol=atol)

        if not sol.success:
            raise RuntimeError("ODE integration failed: " + sol.message)

        return sol.y[:n, -1]


    # def particle_motion_vectorized(self, t, y_flat, params):
    #     """
    #     Vectorized ODE right-hand side.
        
    #     Parameters:
    #     - t: time (scalar)
    #     - y_flat: flattened state array of length 2*n_particles. It is assumed
    #                 that the first n entries are z for each particle and the next n entries are w.
    #     - params: tuple of parameter arrays:
    #         (Cl, rho_f, rho_s, u_star, r, mu, g)
    #         each of shape (n_particles,)
        
    #     Returns:
    #     - dydt_flat: flattened derivative array of the same length.
    #     """
    #     # Determine the number of particles.
    #     n = y_flat.size // 2
    #     # Reshape the state: first row is z, second is w.
    #     y = y_flat.reshape((2, n))
    #     z = y[0, :]   # particle positions
    #     w = y[1, :]   # particle velocities
        
    #     # Unpack parameter arrays.
    #     Cl, rho_f, rho_s, u_star, r, mu, g = params
        
    #     # Prevent division by zero: if z == 0, substitute a small value.
    #     z_fixed = np.where(z == 0, 1e-8, z)
        
    #     # Compute the force terms:
    #     # Lift force term.
    #     lift_term = (0.5 * 3 * Cl * rho_f / (4 * r * rho_s)) * (u_star**2 * (np.log(2))**2) / (z_fixed**2)
    #     # Gravity/buoyancy term.
    #     gravity_term = g * (1 - rho_f / rho_s)
    #     # Drag force term.
    #     drag_term = (9 * mu / (2 * r**2 * rho_s)) * w
        
    #     # Compute derivatives.
    #     dw_dt = lift_term - gravity_term - drag_term
    #     dz_dt = w  # time derivative of z is the velocity w
        
    #     # Stack derivatives into a 2 x n array and flatten it.
    #     dydt = np.vstack((dz_dt, dw_dt))
    #     return dydt.flatten()

    # def calc_resuspension_depth(self, bottom_stress, resuspending):
    #     """
    #     Compute the resuspension depth (i.e. the final particle z position) using a vectorized
    #     ODE integration for all active particles (where resuspending is True). The integration
    #     is performed in one batch call to solve_ivp.
        
    #     For inactive particles (where resuspending is False), the output is set to NaN.
    #     """
    #     gravity = 9.81

    #     # Extract environmental parameters for active particles.
    #     temp = self.environment['sea_water_temperature'][resuspending]
    #     sal = self.environment['sea_water_salinity'][resuspending]
    #     rho_f = self.sea_water_density(temp, sal)

    #     C_l = self.elements.C_l[resuspending]
    #     rho_s = self.elements.rho_s[resuspending]
    #     bot_stress = bottom_stress[resuspending]
        
    #     u_star = np.sqrt(bot_stress / rho_f)
    #     r = self.elements.grain_diameter[resuspending] / 2
    #     mu = self.elements.viscosity_molecular[resuspending]

    #     # Initial conditions for the ODE.
    #     z0_offset = 1e-8
    #     z0_base = 1e-3  # roughness height
    #     z0_initial = z0_base + z0_offset
    #     n_active = np.count_nonzero(resuspending)
    #     print(n_active, 'resuspending particles')
        
    #     # Build initial condition arrays for active particles.
    #     z0_array = np.full(n_active, z0_initial)
    #     w0_array = np.zeros(n_active)
        
    #     # Arrange the state as a 2 x n_active array and flatten it.
    #     y0 = np.vstack((z0_array, w0_array)).flatten()
        
    #     # Integration settings.
    #     npts = 100
    #     t_end = 120  # integration time in seconds
    #     t_span = (0, t_end)
    #     t_eval = np.linspace(t_span[0], t_span[1], npts)
        
    #     # Pack parameters for all active particles into arrays.
    #     # Each parameter array should have shape (n_active,).
    #     params = (C_l, rho_f, rho_s, u_star, r, mu, np.full(n_active, gravity))
        
    #     # Solve the ODE system for all active particles in one call.
    #     sol = solve_ivp(
    #         fun=lambda t, y: self.particle_motion_vectorized(t, y, params),
    #         t_span=t_span, y0=y0, t_eval=t_eval, method='BDF',
    #         rtol=1e-1, atol=1e-5    )
        
    #     if not sol.success:
    #         raise RuntimeError("ODE integration failed: " + sol.message)
        
    #     # The final state is given by sol.y[:, -1] with shape (2*n_active,).
    #     # Reshape it into a 2 x n_active array.
    #     final_state = sol.y[:, -1].reshape((2, n_active))
    #     z_final = final_state[0, :]  # Extract the final z positions for all active particles.
        
    #     return z_final


    def calc_upward_resuspension_velocity(self, bottom_stress, resuspending):
        # Calculate the upwards velocity to give particle getting resuspended.
        # Equation adapted from https://doi.org/10.1061/JYCEAJ.0004937
        E_0 = self.elements.E_0[resuspending]
        porosity = self.elements.porosity[resuspending]
        rho_s = self.elements.rho_s[resuspending]
        tau_crit = self.elements.tau_crit[resuspending]
        bot_stress = bottom_stress[resuspending]
        
        w = ((E_0 * (1-porosity))/rho_s) * ((bot_stress - tau_crit)/(bot_stress))
        return w

    def plot_property_sedimentdrift(self, prop, filename=None, mean=False, labels=None, days=False, num_per_group=None, legend=True):
        """Basic function to plot time series of any element properties."""
        # This function adds onto the plot_property function that the basemodel already provides.
        # If you set days to true, the x axis labels wont include hours and seconds
        # Legend defaults to true. If you have more than one particle per label group.
        # you need to specify num_per_group and labels.
        import matplotlib.pyplot as plt
        from matplotlib import dates

        if not days:
            hfmt = dates.DateFormatter('%d %b %Y %H:%M')
        else:
            hfmt = dates.DateFormatter('%d %b %Y')
        fig = plt.figure()
        ax = fig.gca()
        ax.xaxis.set_major_formatter(hfmt)
        plt.xticks(rotation='vertical')
        start_time = self.start_time
        # In case start_time is unsupported cftime
        start_time = datetime(start_time.year, start_time.month,
                              start_time.day, start_time.hour,
                              start_time.minute, start_time.second)
        times = [
            start_time + n * self.time_step_output
            for n in range(self.steps_output)
        ]
        data = self.history[prop].T[0:len(times), :]
        if mean is True:  # Taking average over elements
            data = np.mean(data, axis=1)
            plt.plot(times, data)
        else:
            if num_per_group is not None:
                colors = plt.cm.viridis(np.linspace(0, 1, int(math.ceil(len(data[0,:]))/num_per_group)))
            for col in range(len(data[0,:])):
                if labels is not None:
                    group = col // num_per_group  # Determine the group (0, 1, 2, or 3)
                    color = colors[group]
                    if col % num_per_group == 0:  # Only label the first line in each group for the legend
                        plt.plot(times, data[:, col], color=color, label=labels[group])
                    else:
                        plt.plot(times, data[:, col], color=color)
                else:                
                    plt.plot(times, data[:,col], label=f"particle_{col}")

            if legend:
                plt.legend()
        plt.title(prop)
        plt.xlabel('Time  [UTC]')
        try:
            plt.ylabel('%s  [%s]' %
                       (prop, self.elements.variables[prop]['units']))
        except:
            plt.ylabel(prop)
        plt.subplots_adjust(bottom=.3)
        plt.grid()
        if filename is None:
            plt.show()
        else:
            plt.savefig(filename)

