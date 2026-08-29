"""Near-bed Rouse resuspension-height sampling.

Covers `opendrift.models.rouse` and the `SedimentDrift.calc_resuspension_height`
dispatcher for the non-ODE modes ('analytic', 'reference', 'turbulent'), plus the
timestep-consistency cap and the production config defaults.
"""
import datetime

import numpy as np

from opendrift.models import rouse
from opendrift.models.sedimentdrift import SedimentDrift

from .sppm_helpers import fabricate_height_model


def test_rouse_number():
    """P = w_s / (beta * kappa * u_star)."""
    assert np.isclose(rouse.rouse_number(1e-3, 0.01), 0.25)
    assert np.isclose(rouse.rouse_number(1e-3, 0.001), 2.5)


def test_powerlaw_bounded_monotonic():
    """Inverse-CDF draw is bounded by [za, zb], monotonic in U, and hugs the bed as P grows."""
    U = np.linspace(0, 1, 7)
    for P in [0.3, 1.0, 2.5]:
        z = rouse.sample_height_powerlaw(np.full_like(U, P), 0.01, 5.0, U)
        assert np.all(np.diff(z) >= -1e-9)
        assert z.min() >= 0.01 - 1e-9 and z.max() <= 5.0 + 1e-9
    z_lo = rouse.sample_height_powerlaw(np.full(5, 0.3), 0.01, 5.0, np.full(5, 0.5))
    z_hi = rouse.sample_height_powerlaw(np.full(5, 2.5), 0.01, 5.0, np.full(5, 0.5))
    assert np.all(z_hi < z_lo)


def test_suspension_threshold_and_loft():
    """P >= P_crit keeps the element in the bedload layer; low P lofts it."""
    h = np.full(4, 100.0)
    z_bed = rouse.resuspension_height(np.full(4, 8e-4), np.full(4, 1e-3), h, 0.01,
                                      seed_layer=5.0, U=np.full(4, 0.9))
    z_loft = rouse.resuspension_height(np.full(4, 5e-2), np.full(4, 1e-3), h, 0.01,
                                       seed_layer=5.0, U=np.full(4, 0.9))
    assert np.allclose(z_bed, 0.01)
    assert np.all(z_loft > 1.0)


def test_reproducible_under_seed():
    """Draws are governed by the global seeded RNG, so runs are reproducible."""
    np.random.seed(42)
    a = rouse.resuspension_height(np.full(3, 0.02), np.full(3, 1e-3),
                                  np.full(3, 100.), 0.01, 5.0)
    np.random.seed(42)
    b = rouse.resuspension_height(np.full(3, 0.02), np.full(3, 1e-3),
                                  np.full(3, 100.), 0.01, 5.0)
    assert np.allclose(a, b)


def test_dispatcher_modes():
    """Each height mode returns sane values, and 'turbulent' orders heights by w_s."""
    o = fabricate_height_model()
    resusp = np.ones(4, dtype=bool)
    bottom_stress = np.full(4, 0.5)          # Pa -> u* ~ 0.022 m/s

    o.set_config('vertical_mixing:resuspension_height_mode', 'analytic')
    z_an = o.calc_resuspension_height(bottom_stress, resusp)
    assert np.all(np.isfinite(z_an)) and np.all(z_an > 0)
    z_an2 = o.calc_resuspension_height(np.full(4, 2.0), resusp)
    assert np.all(z_an2 > z_an)              # equilibrium height scales with u*

    o.set_config('vertical_mixing:resuspension_height_mode', 'reference')
    assert np.allclose(o.calc_resuspension_height(bottom_stress, resusp), 0.01)

    o.set_config('vertical_mixing:resuspension_height_mode', 'turbulent')
    np.random.seed(1)
    z_turb = o.calc_resuspension_height(bottom_stress, resusp)
    assert np.all(z_turb >= 0.01 - 1e-9)
    assert z_turb.max() <= min(100.0, 0.01 + 5.0) + 1e-6
    # slowest settler (w_s=1e-4, smallest P) sits highest; fastest (3e-3) lowest
    mean_z = np.mean([o.calc_resuspension_height(bottom_stress, resusp)
                      for _ in range(2000)], axis=0)
    assert mean_z[3] > mean_z[0] > mean_z[2]


def test_turbulent_seed_capped_by_timestep():
    """The single-step seed never exceeds kappa*u*dt (timestep consistency)."""
    o = fabricate_height_model()
    o.set_config('vertical_mixing:resuspension_height_mode', 'turbulent')
    resusp = np.ones(4, dtype=bool)
    bottom_stress = np.full(4, 2.0)
    rho_f = float(o.sea_water_density(10.0, 35.0))
    u_star = np.sqrt(2.0 / rho_f)
    for dt in [1800.0, 60.0]:
        o.time_step = datetime.timedelta(seconds=dt)
        cap = 0.01 + min(o.get_config('vertical_mixing:resuspension_seed_layer'),
                         rouse.KAPPA * u_star * dt)
        zmax = max(o.calc_resuspension_height(bottom_stress, resusp).max()
                   for _ in range(500))
        assert zmax <= cap + 1e-6


def test_production_defaults():
    """Production defaults: dynamic per-class threshold plus physical Rouse height."""
    o = SedimentDrift(loglevel=50)
    assert o.get_config('vertical_mixing:resuspension_height_mode') == 'turbulent'
    assert o.get_config('vertical_mixing:tau_crit_mode') == 'auto'
    for key, expected in [('vertical_mixing:resuspension_reference_height', 0.01),
                          ('vertical_mixing:resuspension_seed_layer', 5.0),
                          ('vertical_mixing:rouse_beta', 1.0),
                          ('vertical_mixing:rouse_suspension_number', 2.5)]:
        assert np.isclose(o.get_config(key), expected)
