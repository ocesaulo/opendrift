"""Terminal-velocity closures.

Covers `opendrift.models.settling`: the Stokes low-Reynolds limit shared by all
closures, the BB16 shape correction, and a comparison against the measured solid
mineral grains of Maggi (2013) (tests/test_data/sppm/df_particles.csv).
"""
import os

import numpy as np
import pytest

from opendrift.models import settling

RHO_F = 1025.0
MU = 1.4e-3              # dynamic viscosity [kg/m/s], the model default
NU = MU / RHO_F          # kinematic [m^2/s]


def _rms_log(a, b):
    return float(np.sqrt(np.mean((np.log10(a) - np.log10(b)) ** 2)))


def test_stokes_limit():
    """At Re << 1 the stokes, dietrich and spherical bb16 closures coincide."""
    d = np.array([2e-6, 5e-6, 1e-5])
    ws = settling.settling_velocity('stokes', d, 2650.0, RHO_F, NU)
    wd = settling.settling_velocity('dietrich', d, 2650.0, RHO_F, NU)
    wb = settling.settling_velocity('bb16', d, 2650.0, RHO_F, NU, corey_shape_factor=1.0)
    assert np.allclose(ws, wd, rtol=2e-3)
    assert np.allclose(ws, wb, rtol=2e-3)


def test_settling_increases_with_size_and_density():
    """Terminal velocity grows with grain size and with excess density."""
    d = np.array([5e-6, 2e-5, 1e-4, 5e-4])
    w = settling.settling_velocity('dietrich', d, 2650.0, RHO_F, NU)
    assert np.all(np.diff(w) > 0)
    w_light = settling.settling_velocity('dietrich', 1e-4, 1500.0, RHO_F, NU)
    w_heavy = settling.settling_velocity('dietrich', 1e-4, 2650.0, RHO_F, NU)
    assert w_heavy > w_light


def test_buoyant_particle_does_not_sink():
    """A particle lighter than the fluid gets zero settling speed, not a NaN."""
    w = settling.settling_velocity('dietrich', 1e-4, 900.0, RHO_F, NU)
    assert np.all(np.isfinite(w)) and np.all(np.asarray(w) <= 0.0 + 1e-12)


def test_bb16_shape_correction():
    """A sphere matches dietrich at low Re; less spherical grains settle slower."""
    w_sphere = settling.settling_velocity('bb16', 3e-5, 2650.0, RHO_F, NU,
                                          corey_shape_factor=1.0)
    w_die = settling.settling_velocity('dietrich', 3e-5, 2650.0, RHO_F, NU)
    assert np.allclose(w_sphere, w_die, rtol=1e-2)
    previous = None
    for csf in (1.0, 0.8, 0.6, 0.4):
        w = settling.settling_velocity('bb16', 5e-4, 2650.0, RHO_F, NU,
                                       corey_shape_factor=csf)
        if previous is not None:
            assert w < previous
        previous = w


def test_maggi_fractal_reduces_to_solid_grain():
    """Df = 3 with d0 = d makes a fractal aggregate identical to a solid grain."""
    d = np.array([1e-5, 5e-5])
    w_floc = settling.settling_velocity('maggi', d, 2650.0, RHO_F, NU,
                                        d0=d, fractal_dim=3.0)
    w_solid = settling.settling_velocity('stokes', d, 2650.0, RHO_F, NU)
    assert np.allclose(w_floc, w_solid, rtol=5e-2)
    # a genuine fractal aggregate (Df < 3) is lighter, hence slower
    w_fractal = settling.settling_velocity('maggi', d, 2650.0, RHO_F, NU,
                                           d0=1e-6, fractal_dim=2.0)
    assert np.all(w_fractal < w_solid)


def test_vs_maggi_2013_measurements(test_data):
    """Closures reproduce measured settling of solid mineral grains, shape helping."""
    pd = pytest.importorskip('pandas')
    path = os.path.join(test_data, 'sppm', 'df_particles.csv')
    df = pd.read_csv(path)
    df = df[df['delta'] == 3.0].dropna(
        subset=['L_microns', 'v_mm_s', 'rho_s', 'rho_w_fit'])
    df = df[(df['v_mm_s'] > 0) & (df['L_microns'] > 0)]
    assert len(df) > 50, 'fixture unexpectedly small'

    d = df['L_microns'].to_numpy() * 1e-6
    rho_s = df['rho_s'].to_numpy()
    rho_f = df['rho_w_fit'].to_numpy()
    v_meas = df['v_mm_s'].to_numpy() * 1e-3

    rmse_die = _rms_log(settling.settling_velocity(
        'dietrich', d, rho_s, rho_f, MU / rho_f), v_meas)
    rmse_bb = _rms_log(settling.settling_velocity(
        'bb16', d, rho_s, rho_f, MU / rho_f, corey_shape_factor=0.7), v_meas)

    # Reference values at the time of porting: dietrich 0.219, bb16(csf=0.7) 0.202.
    assert rmse_die < 0.25, rmse_die
    assert rmse_bb < 0.25, rmse_bb
    assert rmse_bb < rmse_die, (rmse_bb, rmse_die)
