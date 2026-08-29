"""Dynamic critical-shear-stress closure.

Covers `opendrift.models.floc_strength` and the `SedimentDrift._critical_shear_stress`
dispatcher across tau_crit modes, including back-compatibility of 'constant' and
the growth of cohesive strength with consolidation.
"""
import numpy as np

from opendrift.models import floc_strength as fs

from .sppm_helpers import fabricate_tau_crit_model


def test_shields_curve():
    """Soulsby-Whitehouse Shields: ~150 um fine sand near 0.16 Pa, monotonic in d."""
    rho_f, nu_dyn = 1027.0, 1.4e-3
    d = np.array([10, 63, 150, 300, 1000]) * 1e-6
    tau = fs.tau_crit_shields(d, 2650.0, rho_f, nu_dyn)
    assert np.all(np.isfinite(tau)) and np.all(tau > 0)
    assert np.all(np.diff(tau) > 0)
    assert 0.12 < float(fs.tau_crit_shields(150e-6, 2650.0, rho_f, nu_dyn)) < 0.20


def test_floc_strength_monotonic():
    """Cohesive fractal strength increases with solids fraction phi."""
    phi = np.array([0.05, 0.1, 0.2, 0.4])
    tau = fs.tau_crit_floc(phi, 2.0, 4e-6, 2000.0, 1027.0, c_str=100.0)
    assert np.all(np.diff(tau) > 0)


def test_solids_fraction_bounds():
    """phi = (d/d0)^(Df-3) stays in (0, 1] and falls as the floc grows."""
    phi = fs.solids_fraction(np.array([2e-6, 1e-5, 1e-4]), 1e-6, 2.0)
    assert np.all(phi > 0) and np.all(phi <= 1.0)
    assert np.all(np.diff(phi) < 0)


def test_consolidation_growth():
    """phi relaxes phi0 -> phi_max over the consolidation timescale."""
    phi = fs.consolidate_phi(0.1, np.array([0., 6, 24, 72]) * 3600.0, 86400.0, 0.4)
    assert phi[0] == 0.1
    assert np.all(np.diff(phi) > 0)
    assert phi[-1] < 0.4
    assert fs.consolidate_phi(0.1, 1e12, 86400.0, 0.4) <= 0.4


def test_modes_and_backcompat():
    """'constant' reproduces the prescribed per-element tau_crit; 'auto' routes by class."""
    o = fabricate_tau_crit_model()
    idxs = np.ones(4, dtype=bool)

    o.set_config('vertical_mixing:tau_crit_mode', 'constant')
    assert np.allclose(o._critical_shear_stress(idxs), [0.09, 0.09, 0.03, 0.03])

    o.set_config('vertical_mixing:tau_crit_mode', 'auto')
    tau = o._critical_shear_stress(idxs)
    assert 0.12 < tau[0] < 0.20 and 0.18 < tau[1] < 0.25   # solids via Shields
    assert tau[2] < 0.1 and tau[3] < 0.2                    # cohesives via floc strength


def test_auto_consolidation_raises_cohesive_only():
    """Ageing the deposit stiffens the cohesive elements and leaves solids untouched."""
    o = fabricate_tau_crit_model()
    o.set_config('vertical_mixing:tau_crit_mode', 'auto')
    idxs = np.ones(4, dtype=bool)
    o.elements.time_since_settled[:] = 0.0
    tau0 = o._critical_shear_stress(idxs).copy()
    o.elements.time_since_settled[:] = 72 * 3600.0
    tau1 = o._critical_shear_stress(idxs)
    assert tau1[2] > tau0[2] and tau1[3] > tau0[3]
    assert np.allclose(tau1[:2], tau0[:2])
