"""Synthetic 1DV column forcing for the sediment-physics validation cases.

Builds an in-memory xarray Dataset and hands it to `reader_netCDF_CF_generic`,
so the validation cases need no external forcing file. A column has a flat bed,
a uniform horizontal current, and a prescribed vertical eddy diffusivity profile.

The horizontal grid is deliberately wide: elements advect with the current and
are deactivated on leaving the first reader's bounds, which silently truncates a
run if the domain is too small for `current speed x duration`.
"""
import datetime

import numpy as np
import xarray as xr

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.sedimentdrift import SedimentDrift

KAPPA = 0.4

#: Fields SedimentDrift requires that a bare column does not supply.
_ZERO_FALLBACKS = (
    'land_binary_mask',
    'sea_surface_height',
    'upward_sea_water_velocity',
    'x_wind', 'y_wind',
    'ocean_mixed_layer_thickness',
    'sea_surface_wave_significant_height',
    'sea_surface_wave_period_at_variance_spectral_density_maximum',
    'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment',
    'sea_surface_wave_stokes_drift_x_velocity',
    'sea_surface_wave_stokes_drift_y_velocity',
)


def parabolic_diffusivity(depth, h, u_star, kappa=KAPPA, floor=1e-6):
    """K(z) = kappa u* zb (1 - zb/h), with zb the height above the bed.

    `depth` is positive downward from the surface. This is the profile whose
    steady state, balanced against a constant settling velocity, is the Rouse
    distribution -- so it is the right diffusivity for that benchmark.
    """
    zb = np.clip(h - np.asarray(depth, dtype=float), 0.0, h)
    return np.maximum(kappa * u_star * zb * (1.0 - zb / h), floor)


def rouse_profile(zb, P, h, a):
    """Rouse (1937) equilibrium concentration, normalised to the level `a`."""
    zb = np.asarray(zb, dtype=float)
    return (((h - zb) / zb) * (a / (h - a))) ** P


def column_dataset(h=20.0, u=0.3, v=0.0, u_star=0.02, diffusivity=None,
                   nz=81, hours=48, dt_hours=1.0, n_horiz=5,
                   lon0=0.0, lat0=60.0, dlon=4.0, dlat=2.0,
                   temperature=10.0, salinity=35.0,
                   Hs=0.0, Tp=0.0):
    """In-memory forcing for a flat-bottomed column of depth `h`.

    `diffusivity` may be an array over the nz depth levels (surface first); by
    default the parabolic profile consistent with `u_star` is used.
    """
    z = -np.linspace(0.0, h, nz)
    depth_pos = -z
    K = (parabolic_diffusivity(depth_pos, h, u_star) if diffusivity is None
         else np.broadcast_to(np.asarray(diffusivity, dtype=float), (nz,)))

    times = [datetime.datetime(2000, 1, 1) + datetime.timedelta(hours=dt_hours * i)
             for i in range(int(hours / dt_hours) + 1)]
    lon = lon0 + np.arange(n_horiz) * dlon
    lat = lat0 + np.arange(n_horiz) * dlat
    shape = (len(times), nz, n_horiz, n_horiz)
    ones = np.ones(shape, dtype='f4')

    def _f(value, attrs):
        """Broadcast a scalar, or a per-time series, over the whole grid."""
        value = np.asarray(value, dtype=float)
        if value.ndim == 0:
            field = ones * value
        elif value.shape == (len(times),):
            field = ones * value[:, None, None, None]
        else:
            raise ValueError('expected a scalar or one value per time step, got %s'
                             % (value.shape,))
        return (('time', 'z', 'lat', 'lon'), field.astype('f4'), attrs)

    data = dict(
        u=_f(u, {'standard_name': 'x_sea_water_velocity', 'units': 'm s-1'}),
        v=_f(v, {'standard_name': 'y_sea_water_velocity', 'units': 'm s-1'}),
        temp=_f(temperature, {'standard_name': 'sea_water_temperature',
                              'units': 'degree_Celsius'}),
        salt=_f(salinity, {'standard_name': 'sea_water_salinity', 'units': '1e-3'}),
        K=(('time', 'z', 'lat', 'lon'),
           np.broadcast_to(K[None, :, None, None], shape).astype('f4'),
           {'standard_name': 'ocean_vertical_diffusivity', 'units': 'm2 s-1'}),
        depth=(('lat', 'lon'), np.full((n_horiz, n_horiz), h, dtype='f4'),
               {'standard_name': 'sea_floor_depth_below_sea_level', 'units': 'm'}),
    )
    if Hs > 0:
        data['Hs'] = _f(Hs, {'standard_name': 'sea_surface_wave_significant_height',
                             'units': 'm'})
        data['Tp'] = _f(Tp, {'standard_name': ('sea_surface_wave_period_at_variance'
                                               '_spectral_density_maximum'),
                             'units': 's'})

    return xr.Dataset(
        data_vars=data,
        coords=dict(
            time=('time', times),
            z=('z', z.astype('f4'),
               {'standard_name': 'depth', 'units': 'm', 'positive': 'up'}),
            lat=('lat', lat.astype('f4'),
                 {'standard_name': 'latitude', 'units': 'degrees_north'}),
            lon=('lon', lon.astype('f4'),
                 {'standard_name': 'longitude', 'units': 'degrees_east'}),
        ))


def column_model(dataset, loglevel=50, mixing_timestep=60, **configs):
    """A SedimentDrift wired to `dataset`, with the column defaults applied.

    Extra keyword arguments are config keys with ':' written as '__', e.g.
    ``vertical_mixing__settling_model='bb16'``.
    """
    reader = reader_netCDF_CF_generic.Reader(dataset, name='sppm_column')
    o = SedimentDrift(loglevel=loglevel)
    o.add_reader(reader)
    o.set_config('general:use_auto_landmask', False)
    o.set_config('drift:vertical_mixing', True)
    o.set_config('drift:vertical_advection', False)
    o.set_config('drift:horizontal_diffusivity', 0)
    o.set_config('vertical_mixing:diffusivitymodel', 'environment')
    o.set_config('vertical_mixing:timestep', mixing_timestep)
    for name in _ZERO_FALLBACKS:
        try:
            o.set_config('environment:fallback:' + name, 0)
        except Exception:      # not every field is a declared config
            pass
    for key, value in configs.items():
        o.set_config(key.replace('__', ':'), value)
    return o, reader


def seed_column(o, reader, number, h, w_s, tau_crit=1e-9, z_min=0.01,
                lon=8.0, lat=64.0, rng_seed=1, **kwargs):
    """Seed `number` elements uniformly through the column."""
    np.random.seed(rng_seed)
    o.seed_elements(
        lon=lon, lat=lat, number=number, time=reader.start_time,
        z=-np.random.uniform(z_min, h, number),
        terminal_velocity=np.full(number, -abs(w_s)),
        tau_crit=np.full(number, tau_crit),
        use_stokes=np.zeros(number),
        **kwargs)


def height_above_bed(o, h):
    """Height above the bed of every element still in the domain."""
    return np.clip(h + o.elements.z, 1e-4, h)


def stress_to_speed(tau, rho=1026.95, c_d=0.0021):
    """Near-bed speed giving a quadratic-drag bed stress of `tau` [Pa].

    The 'legacy' bed-stress scheme reduces to tau = rho c_d |u|^2 when the
    deepest velocity cell is thinner than 1 m, so a column with fine vertical
    spacing lets a prescribed stress history be imposed through the current.
    """
    return np.sqrt(np.asarray(tau, dtype=float) / (rho * c_d))


def stress_events(times_hours, events, background=0.0):
    """Piecewise-constant stress history [Pa].

    `events` is a sequence of (start_hour, end_hour, tau) triples.
    """
    tau = np.full(len(times_hours), float(background))
    t = np.asarray(times_hours, dtype=float)
    for start, end, value in events:
        tau[(t >= start) & (t < end)] = value
    return tau
