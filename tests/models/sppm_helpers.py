"""Shared fabrication helpers for the sediment / BBL / settling physics tests.

These build a `SedimentDrift` instance with its element and environment arrays
filled in directly, so the physics routines can be exercised as pure functions
without a reader, a seeded run, or any external forcing file.
"""
import datetime

import numpy as np

from opendrift.models.sedimentdrift import SedimentDrift, SedimentElement

#: Environment fields read by the near-bed / bed-stress routines.
ENV_DTYPE = [
    ('sea_water_temperature', 'f4'),
    ('sea_water_salinity', 'f4'),
    ('sea_floor_depth_below_sea_level', 'f4'),
    ('sea_surface_wave_significant_height', 'f4'),
    ('sea_surface_wave_period_at_variance_spectral_density_maximum', 'f4'),
    ('sea_surface_wave_stokes_drift_x_velocity', 'f4'),
    ('sea_surface_wave_stokes_drift_y_velocity', 'f4'),
]


def make_environment(n, temperature=10.0, salinity=35.0, depth=100.0,
                     Hs=0.0, Tp=0.0, stokes_x=0.0, stokes_y=0.0):
    """Structured environment array of length `n` with the given uniform values."""
    env = np.zeros(n, dtype=ENV_DTYPE)
    env['sea_water_temperature'] = temperature
    env['sea_water_salinity'] = salinity
    env['sea_floor_depth_below_sea_level'] = depth
    env['sea_surface_wave_significant_height'] = Hs
    env['sea_surface_wave_period_at_variance_spectral_density_maximum'] = Tp
    env['sea_surface_wave_stokes_drift_x_velocity'] = stokes_x
    env['sea_surface_wave_stokes_drift_y_velocity'] = stokes_y
    return env


def fabricate_tau_crit_model():
    """Four elements spanning the non-cohesive / cohesive split, for tau_crit tests."""
    o = SedimentDrift(loglevel=50)
    n = 4
    o.elements = SedimentElement(
        lon=np.zeros(n), lat=np.zeros(n), z=-10 * np.ones(n),
        grain_diameter=np.array([150e-6, 300e-6, 4e-6, 10e-6]),
        rho_s=np.array([2650., 2650., 2000., 2000.]),
        sed_class=np.array([0, 0, 1, 1], dtype=np.uint8),
        fractal_dim=np.full(n, 2.0), phi0=np.full(n, 0.1),
        time_since_settled=np.zeros(n),
        tau_crit=np.array([0.09, 0.09, 0.03, 0.03]))
    o.environment = make_environment(n, depth=10.0)
    return o


def fabricate_height_model(n=4, time_step_seconds=1800):
    """Elements spanning a range of settling velocities, for resuspension-height tests."""
    o = SedimentDrift(loglevel=50)
    o.elements = SedimentElement(
        lon=np.zeros(n), lat=np.zeros(n), z=-100.0 * np.ones(n),
        grain_diameter=np.full(n, 30e-6), rho_s=np.full(n, 2650.0),
        C_l=np.full(n, 0.5),
        terminal_velocity=np.array([-1e-3, -1e-3, -3e-3, -1e-4]))
    o.environment = make_environment(n, depth=100.0)
    # the 'turbulent' height mode caps the seed layer at kappa*u*dt, so a real
    # time step must be present (a live run sets this in run()).
    o.time_step = datetime.timedelta(seconds=time_step_seconds)
    return o


def fabricate_bbl_model(n=100, **env_kw):
    """Bare model plus environment, for driving the bed-stress routines directly."""
    o = SedimentDrift(loglevel=50)
    o.elements = SedimentElement(
        lon=np.zeros(n), lat=np.zeros(n), z=-20.0 * np.ones(n),
        grain_diameter=np.full(n, 31e-6), rho_s=np.full(n, 2650.0),
        tau_crit=np.full(n, 0.02))
    o.environment = make_environment(n, **env_kw)
    return o
