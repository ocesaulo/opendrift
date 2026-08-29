"""Bed shear stress and the BBL-validity guard.

`SedimentDrift._bottom_stress_sg2000` routes each element by how well the bottom
boundary layer is resolved by the forcing grid:

  dz_bot > bbl_max_ref_height          -> quadratic drag, tau = rho * c_d * |u|^2
  dz_bot <= cap, no waves at the bed   -> current-only law of the wall
  dz_bot <= cap, waves reach the bed   -> full Styles & Glenn (2000)

These drive the routine directly with fabricated near-bed environments, so no
forcing file is needed.
"""
import numpy as np

from .sppm_helpers import fabricate_bbl_model

U = 0.1  # near-bed current speed [m/s] used throughout


def _env(n, dz_bot, depth):
    return dict(floor_idx=np.arange(n), n_active=n,
                u=np.full(n, U), v=np.zeros(n), visc=np.full(n, 1e-3),
                dz_bot=np.asarray(dz_bot, dtype=float),
                depth=np.full(n, float(depth)))


def _rho(o):
    return float(np.atleast_1d(o.sea_water_density(
        o.environment['sea_water_temperature'][0],
        o.environment['sea_water_salinity'][0])).ravel()[0])


def test_deep_cell_falls_back_to_quadratic_drag():
    """Where the deepest velocity cell is far above the bed, use MITgcm-style drag."""
    n = 20
    o = fabricate_bbl_model(n, depth=500.0)
    o.set_config('vertical_mixing:bbl_scheme', 'sg2000')
    cap = o.get_config('vertical_mixing:bbl_max_ref_height')
    c_d = o.get_config('vertical_mixing:bottom_drag_coefficient')

    tau = o._bottom_stress_sg2000(_env(n, np.full(n, 5 * cap), 500.0))
    assert np.allclose(tau, _rho(o) * c_d * U ** 2, rtol=1e-3)


def test_resolved_cell_uses_log_law():
    """A resolved BBL with no waves gives a current-only log-law stress."""
    n = 20
    o = fabricate_bbl_model(n, depth=500.0)
    o.set_config('vertical_mixing:bbl_scheme', 'sg2000')
    c_d = o.get_config('vertical_mixing:bottom_drag_coefficient')

    tau = o._bottom_stress_sg2000(_env(n, np.full(n, 2.0), 500.0))
    assert np.all(np.isfinite(tau)) and np.all(tau > 0)
    assert not np.allclose(tau, _rho(o) * c_d * U ** 2)


def test_waves_raise_the_stress():
    """Waves reaching a resolved bed give a larger stress than the current alone."""
    n = 20
    o_cur = fabricate_bbl_model(n, depth=20.0)
    o_cur.set_config('vertical_mixing:bbl_scheme', 'sg2000')
    tau_current = o_cur._bottom_stress_sg2000(_env(n, np.full(n, 2.0), 20.0))

    o_wave = fabricate_bbl_model(n, depth=20.0, Hs=2.0, Tp=12.0, stokes_x=0.2)
    o_wave.set_config('vertical_mixing:bbl_scheme', 'sg2000')
    tau_wave = o_wave._bottom_stress_sg2000(_env(n, np.full(n, 2.0), 20.0))

    assert np.all(np.isfinite(tau_wave))
    assert np.nanmean(tau_wave) > np.nanmean(tau_current)


def test_deep_water_waves_do_not_reach_the_bed():
    """The same sea state over deep water leaves the bed stress current-only."""
    n = 20
    o_cur = fabricate_bbl_model(n, depth=2000.0)
    o_cur.set_config('vertical_mixing:bbl_scheme', 'sg2000')
    tau_current = o_cur._bottom_stress_sg2000(_env(n, np.full(n, 2.0), 2000.0))

    o_wave = fabricate_bbl_model(n, depth=2000.0, Hs=2.0, Tp=12.0, stokes_x=0.2)
    o_wave.set_config('vertical_mixing:bbl_scheme', 'sg2000')
    tau_wave = o_wave._bottom_stress_sg2000(_env(n, np.full(n, 2.0), 2000.0))

    assert np.allclose(tau_wave, tau_current, rtol=1e-6)


def test_legacy_scheme_honours_drag_coefficient_config():
    """The legacy bed-stress branch reads bottom_drag_coefficient, not a literal.

    Regression test: that branch previously hard-coded c_d = 0.0021, so changing
    the config silently had no effect on the legacy scheme.
    """
    n = 10
    env = _env(n, np.full(n, 0.5), 500.0)   # thin cell -> no turbulent contribution

    stresses = []
    for c_d in (0.0021, 0.0210):
        o = fabricate_bbl_model(n, depth=500.0)
        o.set_config('vertical_mixing:bbl_scheme', 'legacy')
        o.set_config('vertical_mixing:bottom_drag_coefficient', c_d)
        o._near_bottom_environment = lambda idxs, _e=env: _e
        tau = o.calc_bottom_stress(np.ones(n, dtype=bool))
        stresses.append(float(np.mean(tau[:n])))

    # tenfold c_d must give a tenfold drag stress (the only term for a thin cell)
    assert stresses[1] > stresses[0]
    assert np.isclose(stresses[1] / stresses[0], 10.0, rtol=1e-6), stresses


def test_pickup_velocity_is_callable():
    """The Ariathurai-Partheniades stub reads bed properties from config.

    Regression test: it previously read elements.E_0 / elements.porosity, which do
    not exist on SedimentElement, so any call raised AttributeError.
    """
    n = 5
    o = fabricate_bbl_model(n, depth=100.0)
    resuspending = np.ones(n, dtype=bool)
    bottom_stress = np.full(n, 0.5)          # Pa, well above the seeded tau_crit

    w = o.calc_upward_resuspension_velocity(bottom_stress, resuspending)
    assert np.all(np.isfinite(w)) and np.all(w > 0)

    # scales linearly with the configured erosion-rate constant
    o.set_config('vertical_mixing:erosion_rate_constant', 1e-4)
    w2 = o.calc_upward_resuspension_velocity(bottom_stress, resuspending)
    assert np.allclose(w2 / w, 1e-4 / 5e-5, rtol=1e-6)

    # vanishes as the bed stress approaches the critical stress
    tau_c = float(o.elements.tau_crit[0])
    w_marginal = o.calc_upward_resuspension_velocity(
        np.full(n, tau_c * (1.0 + 1e-9)), resuspending)
    assert np.all(w_marginal < 1e-6 * w2)
