"""Validation against the ROMS/COAWST sediment cases.

Idealized and semi-idealized cases from

  Warner, J. C., Sherwood, C. R., Signell, R. P., Harris, C. K. & Arango, H. G.
  (2008). Development of a three-dimensional, regional, coupled wave, current,
  and sediment-transport model. Computers & Geosciences 34, 1284-1306.

  Sherwood, C. R., Aretxabaleta, A. L., Harris, C. K., Rinehimer, J. P.,
  Verney, R. & Ferre, B. (2018). Cohesive and mixed sediment in the Regional
  Ocean Modeling System (ROMS v3.6) implemented in COAWST. Geosci. Model Dev.
  11, 1849-1871.

These are Eulerian, multi-class, layered-bed models; we are a Lagrangian particle
model with no bed stratigraphy. The cases are therefore adapted, and the
comparison is on physics that both formulations can express: resuspension event
count and intensity, event termination, and the order in which classes deposit.
Where a quantity cannot be compared, the test says so rather than approximating.

See docs/sppmdrift/M1_resuspension_validation.md.
"""
import numpy as np
import pytest

from opendrift.models import floc_strength as fs

from .sppm_helpers import fabricate_tau_crit_model

# Sherwood et al. (2018) section 3.2.1: the four sediment classes.
CLASS_D = np.array([4e-6, 30e-6, 62.5e-6, 140e-6])       # m
CLASS_TAU_CE = np.array([0.05, 0.05, 0.10, 0.10])        # Pa, as prescribed there
CLASS_WS = np.array([0.1e-3, 0.6e-3, 2.0e-3, 8.0e-3])    # m/s

# Sherwood et al. (2018) Fig. 2 configuration.
TAU_CB = 0.10        # bulk critical shear stress of the bed [Pa]
TAU_B = 0.12         # applied bottom stress [Pa]
FC_FULLY_COHESIVE = 0.20   # mud fraction at which the bed is fully cohesive


def sherwood_tau_ce(Pc, tau_cb, tau_c):
    """Sherwood et al. (2018) Eq. 6: tau_ce = max[Pc tau_cb + (1-Pc) tau_c, tau_c].

    Note two properties that our 'mixed' mode does not share: the cohesive term is
    the *bulk bed* stress, identical for every size class, and the result is
    floored at the particle's own critical stress.
    """
    return np.maximum(Pc * tau_cb + (1.0 - Pc) * tau_c, tau_c)


def _reference_rho_f():
    """Seawater density the model itself derives from the fabricated T/S."""
    o = fabricate_tau_crit_model()
    return float(np.atleast_1d(o.sea_water_density(10.0, 35.0)).ravel()[0])


def _shields(d, rho_s=2650.0, rho_f=None, nu_dyn=1.4e-3):
    if rho_f is None:
        rho_f = _reference_rho_f()
    return fs.tau_crit_shields(d, rho_s, rho_f, nu_dyn)


def _model_mixed_tau(fmud, d, rho_s=2650.0):
    """tau_crit from the model's 'mixed' mode for the given bed mud fraction."""
    o = fabricate_tau_crit_model()
    n = len(d)
    o.elements.grain_diameter = np.asarray(d, dtype=np.float32)
    o.elements.rho_s = np.full(n, rho_s, dtype=np.float32)
    o.elements.fractal_dim = np.full(n, 2.0, dtype=np.float32)
    o.elements.phi0 = np.full(n, 0.1, dtype=np.float32)
    o.elements.time_since_settled = np.zeros(n, dtype=np.float32)
    o.elements.sed_class = np.zeros(n, dtype=np.uint8)
    o.elements.tau_crit = np.zeros(n, dtype=np.float32)
    o.set_config('vertical_mixing:tau_crit_mode', 'mixed')
    o.set_config('vertical_mixing:bed_mud_fraction', fmud)
    return o._critical_shear_stress(np.ones(n, dtype=bool))


# --------------------------------------------------------------------------
# S3 - Sherwood et al. (2018) Fig. 2: the mixed-bed closure. No simulation.
# --------------------------------------------------------------------------

def test_s3_pc_ramp_is_monotonic_and_bounded():
    """The cohesive-behaviour parameter ramps monotonically from 0 to 1."""
    lo, hi = 0.05, 0.35
    fmud = np.linspace(0.0, 0.6, 25)
    Pc = np.clip((fmud - lo) / (hi - lo), 0.0, 1.0)
    assert Pc[0] == 0.0 and Pc[-1] == 1.0
    assert np.all(np.diff(Pc) >= 0)


def test_s3_low_mud_limit_reduces_to_shields():
    """At Pc = 0 the effective threshold is the particle's own Shields stress.

    This is the one limit where our formulation and Sherwood Eq. 6 must agree.
    """
    tau_model = _model_mixed_tau(0.0, CLASS_D)
    tau_shields = _shields(CLASS_D)
    assert np.allclose(tau_model, tau_shields, rtol=1e-5)
    assert np.allclose(sherwood_tau_ce(0.0, TAU_CB, tau_shields), tau_shields)


def test_s3_thresholds_increase_with_mud_fraction():
    """Adding mud to the bed never makes a class easier to erode."""
    prev = None
    for fmud in (0.0, 0.1, 0.2, 0.35, 0.5):
        tau = _model_mixed_tau(fmud, CLASS_D)
        if prev is not None:
            assert np.all(tau >= prev - 1e-12), (fmud, tau, prev)
        prev = tau


@pytest.mark.xfail(reason=(
    "Known divergence from Sherwood et al. (2018) Eq. 6, to be resolved in M2. "
    "Our 'mixed' mode blends the particle's own fractal-floc strength, which "
    "scales with the particle's own diameter, instead of the bulk bed stress, "
    "and omits the max(..., tau_c) floor. See docs/sppmdrift/validation/REPORT.md."),
    strict=True)
def test_s3_matches_sherwood_eq6_in_cohesive_limit():
    """In a fully cohesive bed, tau_ce should be max(tau_c, tau_cb) for every class."""
    tau_model = _model_mixed_tau(0.5, CLASS_D)          # Pc = 1
    tau_ref = sherwood_tau_ce(1.0, TAU_CB, _shields(CLASS_D))
    assert np.allclose(tau_model, tau_ref, rtol=0.25), (tau_model, tau_ref)


def test_s3_divergence_is_characterised():
    """Pin the size of the divergence so a fix cannot land silently.

    In the fully cohesive limit our thresholds are far too high for coarse
    classes (the floc-strength law scales with the grain's own diameter) and too
    low for the finest class (no floor at the bulk bed stress). This test records
    the measured behaviour; it must be updated together with any fix.
    """
    tau_model = _model_mixed_tau(0.5, CLASS_D)
    tau_ref = sherwood_tau_ce(1.0, TAU_CB, _shields(CLASS_D))
    ratio = tau_model / tau_ref
    # finest class is too easily eroded; coarse classes far too hard
    assert ratio[0] < 1.0, ratio
    assert ratio[3] > 10.0, ratio
    assert np.all(np.diff(ratio) > 0), ratio      # monotonic in grain size


def test_s3_erosion_ordering_at_applied_stress():
    """Which classes move at the Fig. 2 stress, under the reference formulation.

    Sherwood et al. (2018) Fig. 2 caption states that tau_b = 0.12 Pa is "greater
    than tau_c for clay and silt primary particles, but less than tau_c for sand".
    Our Shields closure must agree: the 140 um class is immobile at this stress
    whatever the mud fraction, and the 62.5 um class sits just below it.
    """
    tau_c = _shields(CLASS_D)
    mobile_noncohesive = TAU_B > sherwood_tau_ce(0.0, TAU_CB, tau_c)
    mobile_cohesive = TAU_B > sherwood_tau_ce(1.0, TAU_CB, tau_c)

    assert mobile_noncohesive[0] and mobile_noncohesive[1]   # clay, silt: move
    assert not mobile_noncohesive[3], tau_c                  # sand: does not
    assert tau_c[3] > TAU_B > TAU_CB

    # made fully cohesive, the fines are governed by the bulk bed stress and are
    # still mobile at 0.12 Pa, while sand remains held by its own Shields stress
    assert mobile_cohesive[0] and mobile_cohesive[1]
    assert not mobile_cohesive[3]


# --------------------------------------------------------------------------
# S1 - Sherwood et al. (2018) section 3.2.1: the double resuspension experiment.
#
# 20 m water depth, one-dimensional in the vertical, flat bottom, no floc
# dynamics. Two bottom-stress events about 1.5 days apart lasting 1.5 and 1 days,
# the first peaking near 1 Pa. Four sediment classes, initially all in the bed.
#
# The published bed is 41 x 1 mm layers holding 25 % of each class; we have no
# bed stratigraphy, so the stratigraphic result (fining-upward grading, a coarse
# lag layer) is compared through the order in which classes leave suspension.
# --------------------------------------------------------------------------

S1_DEPTH = 20.0
S1_HOURS = 120.0
S1_EVENTS = ((6.0, 42.0, 1.00), (78.0, 102.0, 0.45))   # start h, end h, peak Pa


def _s1_stress(hours, dt=1.0):
    """Two half-sine bottom-stress events, as in Fig. 5b."""
    t = np.arange(0.0, hours + dt, dt)
    tau = np.zeros_like(t)
    for t0, t1, peak in S1_EVENTS:
        m = (t >= t0) & (t <= t1)
        tau[m] = peak * np.sin(np.pi * (t[m] - t0) / (t1 - t0))
    return t, tau


@pytest.fixture(scope='module')
def s1_run(tmp_path_factory):
    """Run the double-resuspension case once and share it across the assertions."""
    import datetime
    import xarray as xr
    from .sppm_harness import column_dataset, column_model, stress_to_speed

    n_per_class = 500
    t_h, tau_h = _s1_stress(S1_HOURS)
    ds = column_dataset(h=S1_DEPTH, u=stress_to_speed(tau_h), u_star=0.02,
                        hours=S1_HOURS, dt_hours=1.0, nz=81)
    o, reader = column_model(
        ds,
        vertical_mixing__settling_model='prescribed',
        vertical_mixing__tau_crit_mode='constant',
        vertical_mixing__resuspension_height_mode='turbulent',
        vertical_mixing__bbl_scheme='legacy',
        vertical_mixing__resuspension_seed_layer=5.0)

    np.random.seed(3)
    n = n_per_class * len(CLASS_D)
    o.seed_elements(
        lon=8.0, lat=64.0, number=n, time=reader.start_time, z='seafloor',
        terminal_velocity=np.repeat(-CLASS_WS, n_per_class),
        tau_crit=np.repeat(CLASS_TAU_CE, n_per_class),
        grain_diameter=np.repeat(CLASS_D, n_per_class),
        rho_s=np.full(n, 2650.0), use_stokes=np.zeros(n),
        origin_marker=np.repeat(np.arange(len(CLASS_D)), n_per_class),
        settled=np.ones(n), moving=np.zeros(n))

    outfile = str(tmp_path_factory.mktemp('sppm') / 's1.nc')
    o.run(time_step=300, time_step_output=1800,
          duration=datetime.timedelta(hours=S1_HOURS), outfile=outfile,
          export_variables=['z', 'settled', 'moving', 'origin_marker',
                            'times_resuspended'],
          stop_on_error=True)

    with xr.open_dataset(outfile) as d:
        settled = d['settled'].values
        marker = d['origin_marker'].values[:, 0].astype(int)
        resusp = d['times_resuspended'].values
    hours = np.arange(settled.shape[1]) * 0.5
    return dict(settled=settled, marker=marker, resusp=resusp, hours=hours,
                n_per_class=n_per_class, tau=(t_h, tau_h))


def _suspended_fraction(run, k):
    m = run['marker'] == k
    return np.array([np.nanmean(run['settled'][m, i] == 0)
                     for i in range(run['settled'].shape[1])])


@pytest.mark.slow
def test_s1_no_material_leaves_the_column(s1_run):
    """The domain is wide enough that the experiment conserves particles."""
    lost = np.sum(~np.isfinite(s1_run['settled'][:, -1]))
    assert lost == 0, lost


@pytest.mark.slow
def test_s1_suspended_load_is_ordered_by_grain_size(s1_run):
    """Fines dominate the suspended load; only part of the coarsest class moves.

    Sherwood et al. (2018): "The finer fractions dominated the suspended sediment
    in the water column, which contained only a small fraction of the coarsest
    sand."
    """
    peaks = []
    for k in range(len(CLASS_D)):
        frac = _suspended_fraction(s1_run, k)
        window = (s1_run['hours'] >= 6) & (s1_run['hours'] <= 42)
        peaks.append(frac[window].max())
    peaks = np.array(peaks)
    assert np.all(np.diff(peaks) <= 1e-9), peaks      # non-increasing with size
    assert peaks[0] > 0.95 and peaks[-1] < 0.80, peaks


@pytest.mark.slow
def test_s1_deposition_order_follows_settling_velocity(s1_run):
    """After the stress subsides, classes leave suspension fastest-settling first.

    This is the Lagrangian expression of the fining-upward storm layer and the
    coarse lag: the deposition sequence that builds the graded bed.
    """
    hours = s1_run['hours']
    settle_hour = []
    for k in range(len(CLASS_D)):
        frac = _suspended_fraction(s1_run, k)
        after = np.where((hours >= 42) & (frac < 0.10))[0]
        settle_hour.append(hours[after[0]] if after.size else np.inf)
    # 140 um first, then 62.5, then 30; 4 um never clears within the run
    assert settle_hour[3] < settle_hour[2] < settle_hour[1] < settle_hour[0]
    assert np.isinf(settle_hour[0]), settle_hour


@pytest.mark.slow
def test_s1_finest_class_still_suspended_after_five_days(s1_run):
    """Sherwood: "The finest material (4 um) remained mostly in suspension"."""
    frac = _suspended_fraction(s1_run, 0)
    assert frac[-1] > 0.5, frac[-1]
    # and the two sand classes are fully deposited by then
    assert _suspended_fraction(s1_run, 2)[-1] < 0.05
    assert _suspended_fraction(s1_run, 3)[-1] < 0.05


@pytest.mark.slow
def test_s1_second_event_does_not_suppress_coarse_resuspension(s1_run):
    """DOCUMENTED FAILURE: a weaker second event should move much less sand.

    Sherwood et al. (2018): "The second stress pulse eroded the bed down to 1 cm
    but only resuspended minimal amounts of the 140 um sand." Their erosion is a
    flux, E = E_0 (1 - phi) (tau_b/tau_ce - 1), so a weaker event moves
    proportionally less material.

    Ours is an all-or-nothing threshold crossing: once tau_b exceeds tau_ce every
    settled element of that class lifts off, however small the excess. At the
    second event's 0.45 Pa peak against a 0.10 Pa threshold, essentially all the
    sand is resuspended again.

    This test pins the defect rather than hiding it. Resolving it means adopting
    a probabilistic pickup flux (M2/B4); when that lands, this test should be
    replaced by the suppression assertion in its docstring.
    """
    frac = _suspended_fraction(s1_run, 3)          # 140 um sand
    hours = s1_run['hours']
    peak1 = frac[(hours >= 6) & (hours <= 42)].max()
    peak2 = frac[(hours >= 78) & (hours <= 102)].max()
    # what Sherwood would give, and what we should assert once the flux exists:
    #     assert peak2 < 0.25 * peak1
    assert peak2 > 0.75 * peak1, (peak1, peak2)


@pytest.mark.slow
def test_s1_coarse_classes_hop_repeatedly(s1_run):
    """DOCUMENTED FAILURE: settled elements above threshold re-lift every step.

    A settled element whose critical stress is exceeded is lifted again on every
    time step, so the coarse classes accumulate hundreds of resuspension events
    within two stress events. This is the mechanism behind the ~45 % out-of-domain
    loss of the fastest settling classes reported in the parent DDT study; here it
    is isolated in a case with a known answer.

    Two events should give at most a few lift-offs per particle.

    Note: this diagnostic was itself unreliable until `times_resuspended` was
    widened from uint8, which silently wrapped at 255.
    """
    per_particle = []
    for k in range(len(CLASS_D)):
        m = s1_run['marker'] == k
        last = np.nanmax(np.where(np.isfinite(s1_run['resusp'][m]),
                                  s1_run['resusp'][m], np.nan), axis=1)
        per_particle.append(np.nanmean(last))
    per_particle = np.array(per_particle)
    assert np.all(np.diff(per_particle) > 0), per_particle   # worse for coarser
    assert per_particle[3] > 100.0, per_particle             # ~299 measured


@pytest.mark.slow
def test_s1_lift_off_rate_does_not_converge():
    """DOCUMENTED FAILURE, and M2's acceptance gate.

    The decisive property. A physical erosion rate must converge under timestep
    refinement: halving dt must not change how much material leaves the bed.
    Ours does change, because `resuspension()` lifts every settled element above
    threshold once per outer step, so the lift-off count is set by how often the
    routine is called rather than by the flow.

    Measured over a single stress event on one class, varying only dt:

        dt = 900 s -> 67.5 lift-offs;  450 s -> 91.7;  225 s -> 110.3

    a +63 % change over a 4x refinement, monotonic, with no sign of a limit.
    Since each lift-off lofts a grain into faster water before it settles back,
    the horizontal transport inherits the same dependence.

    This test asserts the *defect*. When the pickup flux of M2/B4 lands it must
    start failing, and should then be inverted to assert convergence (the eroded
    mass agreeing across dt to within sampling noise). Do not delete it.
    """
    import datetime
    from .sppm_harness import column_dataset, column_model, stress_to_speed

    depth, hours, n = 20.0, 30.0, 200
    t = np.arange(0, hours + 1, 1.0)
    tau = np.where((t >= 3) & (t <= 27), np.sin(np.pi * (t - 3) / 24), 0.0)

    counts = {}
    for dt in (900, 225):
        ds = column_dataset(h=depth, u=stress_to_speed(tau), u_star=0.02,
                            hours=hours, dt_hours=1.0, nz=81)
        o, reader = column_model(
            ds,
            vertical_mixing__settling_model='prescribed',
            vertical_mixing__tau_crit_mode='constant',
            vertical_mixing__resuspension_height_mode='turbulent',
            vertical_mixing__bbl_scheme='legacy',
            vertical_mixing__resuspension_seed_layer=5.0)
        np.random.seed(3)
        o.seed_elements(
            lon=8.0, lat=64.0, number=n, time=reader.start_time, z='seafloor',
            terminal_velocity=np.full(n, -8.0e-3), tau_crit=np.full(n, 0.10),
            grain_diameter=np.full(n, 140e-6), rho_s=np.full(n, 2650.0),
            use_stokes=np.zeros(n), settled=np.ones(n), moving=np.zeros(n))
        o.run(time_step=dt, time_step_output=3600,
              duration=datetime.timedelta(hours=hours), stop_on_error=True)
        counts[dt] = float(np.mean(o.elements.times_resuspended))

    # a convergent scheme would agree within sampling noise; this one grows
    growth = counts[225] / counts[900]
    assert growth > 1.3, counts
