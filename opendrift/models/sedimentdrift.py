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
from joblib import Parallel, delayed
from multiprocessing import cpu_count
import math


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


        # By default, sediments do not strand towards coastline
        # TODO: A more sophisticated stranding algorithm is needed
        self._set_config_default('general:coastline_action', 'previous')

        # Vertical mixing is enabled as default
        self._set_config_default('drift:vertical_mixing', True)

        self.use_parallel = True

    def update(self):
        """Update positions and properties of sediment particles.
        """
        # Advecting here all elements, but want to soon add
        # possibility of not moving settled elements, until
        # they are resuspended. May then need to send a boolean
        # array to advection methods below
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

        self.resuspension()
        self.deactivate_elements(self.elements.beached == 1, reason='beached')
        self.remove_deactivated_elements()

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
            #  self.get_distance_cell_center_above(settling)  # no longer needed here

    def resuspension(self):
        """Resuspending elements when critical shear stress is exceed"""
        ## Old resuspension condition
        # threshold = self.get_config('vertical_mixing:resuspension_threshold')
        # resuspending = np.logical_and(self.current_speed() > threshold, self.elements.moving==0)

        threshold = self.elements.tau_crit  # tau_crit should become a function of the element and environment

        settled = np.logical_and(self.elements.moving == 0, self.elements.settled == 1)
        if np.sum(settled) > 0:
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
            # Calculate the particle z position from height above bottom following particle motion scheme.
            height_above_bottom = self.calc_resuspension_depth(bottom_stress, resuspending)
            # print(np.max(height_above_bottom))
            new_z = self.elements.z[resuspending] + height_above_bottom
            new_z[new_z > 0] = 0
            # Update the z position of the particle
            self.elements.z[resuspending] = new_z
            # Switch particles that are resuspending out of settled mode
            self.elements.settled[resuspending] = 0
            # keep track of the latest resuspension height:
            self.elements.latest_resuspension_height[resuspending] = height_above_bottom

    def calc_bottom_stress(self, idxs):
        """
        Compute bottom stress for a set of particles (or grid cells) identified by the boolean mask `idxs`.
        The returned array has the same shape as the full elements array. For indices not included in `idxs`,
        the bottom stress is set to NaN.
        
        This function uses a vectorized search (via np.searchsorted) to locate the appropriate vertical
        cell for each active particle, then computes the bottom stress with the same logic as before.
        """
        # 1. Extract data only for active particles (where idxs is True)
        floor_idx = np.where(idxs)[0]  # indices in the full array that are active
        if floor_idx.size == 0:
            # No active particles: return an array of zeros with the full shape.
            return np.zeros(self.elements.z.shape)
            
        zs   = self.elements.z[idxs]
        lons = self.elements.lon[idxs]
        lats = self.elements.lat[idxs]
        atimes = self.time

        variables = ['x_sea_water_velocity', 'y_sea_water_velocity',
                    'ocean_vertical_diffusivity', 'sea_floor_depth_below_sea_level']
        profiles = ['x_sea_water_velocity', 'y_sea_water_velocity', 'ocean_vertical_diffusivity']

        # Get environment data (both cell-centered and profile)
        env_var, env_profs, amiss = self.env.get_environment(variables, atimes, lons, lats, zs, profiles)

        # seafloor depth (convert sign so that depth is positive downward)
        # (Note: self.environment['sea_floor_depth_below_sea_level'] is assumed to be defined over the entire domain)
        el_seafloor_depth_z = -1 * self.environment['sea_floor_depth_below_sea_level'][idxs]

        # Get the vertical grid from the profile data.
        # Here we assume that zgrid is a 1D array (e.g., shape (n_levels,)) that is monotonic.
        zgrid = env_profs['z']  # e.g., descending (surface to seafloor)

        # 2. Initialize arrays for active particles
        n_active = floor_idx.size
        # For each active particle, we want to determine the local near-bottom cell.
        # Set minimum bottom layer thickness and a critical layer thickness.
        dz_bot_min = 0.1                # minimum allowed thickness (m)
        crit_last_layer_thickness = 10.  # threshold for turbulent contribution

        # Allocate arrays of length n_active for near-bottom quantities.
        # (They will be filled with vectorized operations below.)
        u_bot    = np.empty(n_active)
        v_bot    = np.empty(n_active)
        visc_bot = np.empty(n_active)
        dz_bot   = np.full(n_active, dz_bot_min)
        turb_contrib = np.zeros(n_active)

        # 3. Find the near-bottom level for each active particle
        #
        # For each active particle with seafloor depth `adep` (from el_seafloor_depth_z),
        # we need the index of the last grid cell in zgrid that is above the seafloor.
        #
        # The original code did:
        #     last_indz = np.argwhere(zgrid > adep)
        # and used last_indz[-1] if any exist; otherwise, fallback to index 0.
        #
        # We can do this vectorially if we assume that zgrid is monotonic. However,
        # since np.searchsorted requires an ascending array, we invert the sign.
        #
        # For each active particle, we compute:
        #    search_idx = np.searchsorted(-zgrid, -adep, side='right')
        # If search_idx > 0 then the desired index is search_idx - 1; otherwise, we use 0.
        search_indices = np.searchsorted(-zgrid, -el_seafloor_depth_z, side='right')
        selected_idx = np.where(search_indices > 0, search_indices - 1, 0)

        # 4. Extract near-bottom values from the profile arrays.
        # We assume that env_profs['x_sea_water_velocity'], etc., have shape (n_levels, n_active)
        particle_range = np.arange(n_active)
        u_bot    = env_profs['x_sea_water_velocity'][selected_idx, particle_range]
        v_bot    = env_profs['y_sea_water_velocity'][selected_idx, particle_range]
        visc_bot = env_profs['ocean_vertical_diffusivity'][selected_idx, particle_range]

        # Compute the local layer thickness for each active particle.
        dz_bot = zgrid[selected_idx] - el_seafloor_depth_z

        # Log a debug message for any particle that got the fallback index.
        fallback = (search_indices == 0)
        if np.any(fallback):
            for i in np.where(fallback)[0]:
                self.elements.beached[floor_idx[i]] = 1
                logger.debug('Element %s settled at seafloor shallower than shallowest grid cell' %
                            self.elements.ID[floor_idx[i]])

        # 5. Compute the bottom stress using the same formulas as before
        r_b = 0
        c_d = 0.0021  # drag coefficient (to be updated later for a log-wall model)
        _2KE = u_bot**2 + v_bot**2  # twice the kinetic energy (ignoring vertical velocity)
        vel_mag = np.sqrt(_2KE)

        # Turbulent contribution is computed only if the layer thickness is above a threshold.
        turb_layers = dz_bot >= crit_last_layer_thickness
        turb_contrib[turb_layers] = 2 * visc_bot[turb_layers] / dz_bot[turb_layers] * vel_mag[turb_layers]

        # Mean-flow contribution (as in the original code)
        mean_flow_contrib = (r_b + c_d * np.sqrt(_2KE)) * vel_mag
        bottom_stress_pseudo_energy_magnitude = mean_flow_contrib + turb_contrib

        # Multiply by the sea water density.
        computed_bottom_stress = self.sea_water_density() * bottom_stress_pseudo_energy_magnitude

        # 6. Place the computed values back into an array of full size
        full_bottom_stress = np.zeros(self.elements.z.shape)
        # For the indices in the full array corresponding to active particles, insert the computed values.
        full_bottom_stress[floor_idx] = computed_bottom_stress

        return full_bottom_stress

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


    def get_optimal_n_jobs(fraction=0.75, max_jobs=None):
        available = cpu_count()
        jobs = int(available * fraction)
        return min(jobs, max_jobs) if max_jobs else jobs

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


    def solve_batch(batch, grain_diameter, bottom_stress, rtol, atol, dameth):
        gravity = 9.81
        C_l = 0.5
        rho_f = 1026.
        rho_s = 2650.
        mu = 1.4e-3

        r = grain_diameter[batch] / 2
        u_star = np.sqrt(bottom_stress[batch] / rho_f)

        C_l_b = np.full_like(r, C_l)
        rho_f_b = np.full_like(r, rho_f)
        rho_s_b = np.full_like(r, rho_s)
        mu_b = np.full_like(r, mu)
        g_b = np.full_like(r, gravity)

        z0 = np.full(len(batch), 1e-3 + 1e-8)
        w0 = np.zeros(len(batch))
        y0 = np.vstack((z0, w0)).flatten()
        params = (C_l_b, rho_f_b, rho_s_b, u_star, r, mu_b, g_b)

        sol = solve_ivp(
            fun=lambda t, y: particle_motion_vectorized(t, y, params),
            t_span=(0, 120), y0=y0, method=dameth,
            rtol=rtol, atol=atol
        )

        if not sol.success:
            raise RuntimeError("ODE integration failed: " + sol.message)

        final_state = sol.y[:, -1].reshape((2, len(batch)))
        return batch, final_state[0, :]

    def calc_resuspension_depth(grain_diameter, bottom_stress, n_resuspending, rtol=1e-6, atol=1e-8, dameth='BDF', batch_size=500, n_jobs=-1):
        if n_resuspending == 0:
            return np.array([])

        # Ensure inputs are arrays of length n_resuspending
        grain_diameter = np.broadcast_to(grain_diameter, (n_resuspending,)) if np.isscalar(grain_diameter) else np.asarray(grain_diameter)
        bottom_stress = np.broadcast_to(bottom_stress, (n_resuspending,)) if np.isscalar(bottom_stress) else np.asarray(bottom_stress)

        active_indices = np.arange(n_resuspending)
        batches = [active_indices[i:i+batch_size] for i in range(0, n_resuspending, batch_size)]

        results = Parallel(n_jobs=n_jobs, prefer='processes')(
            delayed(solve_batch)(batch, grain_diameter, bottom_stress, rtol, atol, dameth) for batch in batches
        )

        z_final_all = np.full(n_resuspending, np.nan)
        for idx_local, z_local in results:
            z_final_all[idx_local] = z_local

        return z_final_all

    # def particle_motion(self, t, y, Cl, rho_f, rho_s, u_star, r, mu, g):
    #     z, w = y
    #     if z == 0e0:
    #         z = 1e-8  # Prevent division by zero

    #     # Lift force
    #     lift_term = (.5 * 3 * Cl * rho_f / (4 * r * rho_s)) * (u_star**2 * (np.log(2))**2) / z**2
    #     # Gravity/Buoyancy force
    #     gravity_term = g * (1 - rho_f / rho_s)
    #     # Drag force
    #     drag_term = (9 * mu / (2 * r**2 * rho_s)) * w
    #     # Acceleration
    #     dw_dt = lift_term - gravity_term - drag_term
    #     dz_dt = w
    #     return np.array([dz_dt, dw_dt])


    # def solve_particle_motion(self, params, y0, t_span, t_eval):
    #     """ Wrapper to call solve_ivp for each particle. """
    #     solution = solve_ivp(self.particle_motion, t_span, y0, t_eval=t_eval,
    #                         args=params, method="BDF", rtol=1e-6, atol=1e-8)
    #     return solution.y[0][-1], solution.y[1][-1]


    # def calc_resuspension_depth(self, bottom_stress, resuspending):
    #     gravity = 9.81

    #     temp = self.environment['sea_water_temperature'][resuspending]
    #     sal = self.environment['sea_water_salinity'][resuspending]
    #     rho_f = self.sea_water_density(temp, sal)

    #     C_l = self.elements.C_l[resuspending]
    #     rho_s = self.elements.rho_s[resuspending]
    #     bot_stress = bottom_stress[resuspending]
        
    #     u_star = np.sqrt(bot_stress / rho_f)
    #     r = self.elements.grain_diameter[resuspending] / 2
    #     mu = self.elements.viscosity_molecular[resuspending]

    #     z0 = 1e-3  # Roughness height
    #     z0_initial = z0 + 1e-8
    #     w0_initial = 0
    #     y0 = [z0_initial, w0_initial]

    #     npts = 100
    #     t_end = 120
    #     t_span = (0, t_end)
    #     t_eval = np.linspace(t_span[0], t_span[1], npts)

    #     params_list = [(C_l[i], rho_f[i], rho_s[i], u_star[i], r[i], mu[i], gravity) for i in range(len(rho_f))]

    #     if self.use_parallel:
    #         results = Parallel(n_jobs=-1, prefer="threads")(
    #             delayed(self.solve_particle_motion)(params, y0, t_span, t_eval) for params in params_list
    #         )
    #     else:
    #         results = [self.solve_particle_motion(params, y0, t_span, t_eval) for params in params_list]

    #     zend, wend = zip(*results)
    #     return np.array(zend)

    def particle_motion_vectorized(self, t, y_flat, params):
        """
        Vectorized ODE right-hand side.
        
        Parameters:
        - t: time (scalar)
        - y_flat: flattened state array of length 2*n_particles. It is assumed
                    that the first n entries are z for each particle and the next n entries are w.
        - params: tuple of parameter arrays:
            (Cl, rho_f, rho_s, u_star, r, mu, g)
            each of shape (n_particles,)
        
        Returns:
        - dydt_flat: flattened derivative array of the same length.
        """
        # Determine the number of particles.
        n = y_flat.size // 2
        # Reshape the state: first row is z, second is w.
        y = y_flat.reshape((2, n))
        z = y[0, :]   # particle positions
        w = y[1, :]   # particle velocities
        
        # Unpack parameter arrays.
        Cl, rho_f, rho_s, u_star, r, mu, g = params
        
        # Prevent division by zero: if z == 0, substitute a small value.
        z_fixed = np.where(z == 0, 1e-8, z)
        
        # Compute the force terms:
        # Lift force term.
        lift_term = (0.5 * 3 * Cl * rho_f / (4 * r * rho_s)) * (u_star**2 * (np.log(2))**2) / (z_fixed**2)
        # Gravity/buoyancy term.
        gravity_term = g * (1 - rho_f / rho_s)
        # Drag force term.
        drag_term = (9 * mu / (2 * r**2 * rho_s)) * w
        
        # Compute derivatives.
        dw_dt = lift_term - gravity_term - drag_term
        dz_dt = w  # time derivative of z is the velocity w
        
        # Stack derivatives into a 2 x n array and flatten it.
        dydt = np.vstack((dz_dt, dw_dt))
        return dydt.flatten()

    def calc_resuspension_depth(self, bottom_stress, resuspending):
        """
        Compute the resuspension depth (i.e. the final particle z position) using a vectorized
        ODE integration for all active particles (where resuspending is True). The integration
        is performed in one batch call to solve_ivp.
        
        For inactive particles (where resuspending is False), the output is set to NaN.
        """
        gravity = 9.81

        # Extract environmental parameters for active particles.
        temp = self.environment['sea_water_temperature'][resuspending]
        sal = self.environment['sea_water_salinity'][resuspending]
        rho_f = self.sea_water_density(temp, sal)

        C_l = self.elements.C_l[resuspending]
        rho_s = self.elements.rho_s[resuspending]
        bot_stress = bottom_stress[resuspending]
        
        u_star = np.sqrt(bot_stress / rho_f)
        r = self.elements.grain_diameter[resuspending] / 2
        mu = self.elements.viscosity_molecular[resuspending]

        # Initial conditions for the ODE.
        z0_offset = 1e-8
        z0_base = 1e-3  # roughness height
        z0_initial = z0_base + z0_offset
        n_active = np.count_nonzero(resuspending)
        
        # Build initial condition arrays for active particles.
        z0_array = np.full(n_active, z0_initial)
        w0_array = np.zeros(n_active)
        
        # Arrange the state as a 2 x n_active array and flatten it.
        y0 = np.vstack((z0_array, w0_array)).flatten()
        
        # Integration settings.
        npts = 100
        t_end = 120  # integration time in seconds
        t_span = (0, t_end)
        t_eval = np.linspace(t_span[0], t_span[1], npts)
        
        # Pack parameters for all active particles into arrays.
        # Each parameter array should have shape (n_active,).
        params = (C_l, rho_f, rho_s, u_star, r, mu, np.full(n_active, gravity))
        
        # Solve the ODE system for all active particles in one call.
        sol = solve_ivp(
            fun=lambda t, y: self.particle_motion_vectorized(t, y, params),
            t_span=t_span, y0=y0, t_eval=t_eval, method='BDF',
            rtol=1e-6, atol=1e-8    )
        
        if not sol.success:
            raise RuntimeError("ODE integration failed: " + sol.message)
        
        # The final state is given by sol.y[:, -1] with shape (2*n_active,).
        # Reshape it into a 2 x n_active array.
        final_state = sol.y[:, -1].reshape((2, n_active))
        z_final = final_state[0, :]  # Extract the final z positions for all active particles.
        
        return z_final

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

