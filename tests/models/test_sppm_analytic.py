"""Analytic benchmarks for the settling / vertical-mixing / resuspension chain.

These have closed-form answers, need no external forcing, and isolate one
process each. See docs/sppmdrift/M1_resuspension_validation.md, Tier A.
"""
import datetime

import numpy as np
import pytest

from .sppm_harness import (KAPPA, column_dataset, column_model, seed_column,
                           height_above_bed, parabolic_diffusivity)

H = 20.0            # column depth [m]
U_STAR = 0.02       # friction velocity used for the parabolic diffusivity [m/s]
ZA = 0.01           # resuspension reference height [m]


def test_parabolic_diffusivity_shape():
    """K vanishes at both boundaries and peaks at mid-depth."""
    depth = np.linspace(0.0, H, 101)
    K = parabolic_diffusivity(depth, H, U_STAR)
    assert K[0] < 1e-5 and K[-1] < 1e-5
    assert np.isclose(depth[np.argmax(K)], H / 2, atol=H / 20)
    assert np.isclose(K.max(), KAPPA * U_STAR * H / 4, rtol=1e-2)


@pytest.mark.slow
def test_settling_column_deposition_rate():
    """With no turbulence, the settled fraction is min(1, w_s t / h).

    A cloud spread uniformly through a column of depth h and settling at a
    constant w_s deposits linearly in time: after t, everything initially within
    w_s*t of the bed has landed.
    """
    w_s = 2.0e-3                      # m/s -> full column swept in h/w_s = 10000 s
    hours = 2.0                       # 7200 s, so ~72 % expected
    ds = column_dataset(h=H, u=0.0, u_star=U_STAR, hours=hours + 1)
    o, reader = column_model(
        ds,
        drift__vertical_mixing=False,           # isolate settling
        vertical_mixing__settling_model='prescribed',
        vertical_mixing__tau_crit_mode='constant')
    # tau_crit huge so nothing is resuspended once it lands
    seed_column(o, reader, 4000, H, w_s, tau_crit=1e6)
    o.run(time_step=60, time_step_output=3600,
          duration=datetime.timedelta(hours=hours), stop_on_error=True)

    expected = min(1.0, w_s * hours * 3600.0 / H)
    got = np.sum(o.elements.settled == 1) / len(o.elements.settled)
    assert abs(got - expected) < 0.05, (got, expected)


@pytest.mark.slow
def test_well_mixed_condition():
    """A neutrally buoyant cloud stays uniform in a strongly non-uniform K(z).

    This is the standard consistency requirement on a vertical random walk
    (Thomson 1987; Visser 1997; Ross & Sharples 2004): with no settling, an
    initially uniform distribution must remain uniform even though K varies by
    orders of magnitude over the column. The classic failure mode is spurious
    accumulation where K is small, i.e. at both boundaries here.

    Only the near-bed tenth is excluded: the bed is not a reflecting boundary in
    this model but a settle-and-resuspend cycle that re-releases elements at a
    fixed reference height, which is a real boundary treatment but not the one
    this condition is posed for. The surface *is* a reflecting boundary, so the
    upper column -- where K collapses by four orders of magnitude -- is fair game
    and is the more demanding of the two checks.
    """
    hours = 12
    ds = column_dataset(h=H, u=0.0, u_star=U_STAR, hours=hours + 1)
    o, reader = column_model(
        ds,
        vertical_mixing__settling_model='prescribed',
        vertical_mixing__tau_crit_mode='constant',
        vertical_mixing__resuspension_height_mode='reference',
        vertical_mixing__resuspension_reference_height=ZA)
    seed_column(o, reader, 6000, H, w_s=0.0, tau_crit=1e-9)
    o.run(time_step=120, time_step_output=3600,
          duration=datetime.timedelta(hours=hours), stop_on_error=True)

    zb = height_above_bed(o, H)

    # Two regions. The interior is the weak check -- excluding the outer tenths
    # also excludes most of the variation in K, which only spans ~3x there. The
    # upper column is the strong one: the surface is a genuine reflecting
    # boundary, so it can be probed right up to it, and K falls by four orders of
    # magnitude across the top 30 % -- exactly where a biased scheme would pile
    # elements up.
    for lo_f, hi_f in [(0.10, 0.90), (0.70, 1.00)]:
        lo, hi = lo_f * H, hi_f * H
        counts, _ = np.histogram(zb[(zb > lo) & (zb < hi)], bins=8, range=(lo, hi))
        expected = counts.sum() / len(counts)
        # Poisson noise on a few hundred per bin is ~5 %; allow 20 % excursions
        assert np.all(np.abs(counts - expected) / expected < 0.20), (lo_f, counts)


@pytest.mark.slow
@pytest.mark.parametrize('P', [0.25, 0.5])
def test_rouse_equilibrium_profile(P):
    """The suspended profile approaches Rouse for a settling-mixing balance.

    With K(z) = kappa u* z (1 - z/h) and a constant settling velocity, the steady
    state of the vertical advection-diffusion balance is the Rouse distribution
    with exponent P = w_s / (kappa u*). Recovering P from the modelled profile
    tests settling, the random-walk mixing and the resuspension release height
    jointly.

    'reference' release is used rather than the production 'turbulent' mode: the
    latter draws lift-off heights from the Rouse distribution itself, so testing
    a Rouse profile against it would be partly circular.

    Accuracy is limited by the settle-and-resuspend cycle. Elements resting on
    the bed awaiting the next resuspension are outside the suspended profile, and
    their number scales with the outer time step, so the recovered exponent is
    biased low. See docs/sppmdrift/validation/REPORT.md.
    """
    w_s = P * KAPPA * U_STAR
    hours = 24
    ds = column_dataset(h=H, u=0.3, u_star=U_STAR, hours=hours + 2)
    o, reader = column_model(
        ds,
        vertical_mixing__settling_model='prescribed',
        vertical_mixing__tau_crit_mode='constant',
        vertical_mixing__resuspension_height_mode='reference',
        vertical_mixing__resuspension_reference_height=ZA,
        vertical_mixing__bbl_scheme='legacy')
    seed_column(o, reader, 4000, H, w_s, tau_crit=1e-9)
    o.run(time_step=300, time_step_output=3600,
          duration=datetime.timedelta(hours=hours), stop_on_error=True)

    zb = height_above_bed(o, H)
    edges = np.logspace(np.log10(ZA), np.log10(H), 25)
    counts, _ = np.histogram(zb, bins=edges)
    centres = np.sqrt(edges[1:] * edges[:-1])
    concentration = counts / np.diff(edges)
    ok = counts > 20
    assert ok.sum() >= 8, counts

    P_fit = np.polyfit(np.log((H - centres[ok]) / centres[ok]),
                       np.log(concentration[ok]), 1)[0]
    # measured 0.223 (P=0.25) and 0.404 (P=0.5): biased low, monotonic in P
    assert 0.55 * P < P_fit < 1.15 * P, (P, P_fit)
