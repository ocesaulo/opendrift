# M1 — Resuspension physics: validation campaign

## 1. Goal & scope

The resuspension physics is **already built**. What is missing is evidence that
it is right. M1's deliverable is therefore a **benchmark suite**, not new
closures.

"Done" means: a runnable, mostly CI-fast suite in which every case states a
metric, a tolerance, a published value and a model value, with honest fails
recorded rather than tuned away; plus a written verdict on whether the
resuspension physics is fit for the production runs that depend on it.

**Out of scope.** New physics. If a case fails because a scheme is wrong, that is
an M1 *finding* and the fix belongs to M2 — with one flagged exception (§7,
pickup flux). Also out of scope: a real wave forcing reader (WW3/ERA5); M1's
wave path uses synthetic forcing, which is sufficient to validate the physics and
does not block on data engineering.

## 2. Current state

Built, and defaulted on in production since 2026-06-19:

| Piece | Where | Default |
|---|---|---|
| Bed stress, wave–current BBL | `sedimentdrift.py:983 _bottom_stress_sg2000` | `bbl_scheme='sg2000'` |
| BBL-validity guard | same, routes on `dz_bot` vs `bbl_max_ref_height` | 10 m |
| Dynamic erosion threshold | `sedimentdrift.py:768 _critical_shear_stress` | `tau_crit_mode='auto'` |
| Resuspension height | `sedimentdrift.py:1272 calc_resuspension_height` | `resuspension_height_mode='turbulent'` |
| Driver | `sedimentdrift.py:838 resuspension` | — |

Validated so far (all now in `tests/`, see M0):

- Styles & Glenn (2000) Table 2, max abs error 0.000 — the *only* published
  benchmark the bed-stress scheme has ever been checked against.
- Settling vs Maggi (2013) measured grains, log10-RMSE 0.20.
- Internal consistency: monotonicity, bounds, reproducibility, config defaults.

**Never checked:** the equilibrium suspended-sediment profile, the well-mixed
condition, any erosion/deposition benchmark, and any cross-model comparison. The
`turbulent` height scheme and the `auto` threshold became production defaults
without a full regression — `RESUSPENSION_MECHANICS_PLAN.md` says so explicitly
("A full bit-for-bit `resusp09` regression via `ddt_sim` remains a good gate
before switching defaults"; the defaults were switched anyway).

## 3. Literature gate — CLEARED

Two papers, both on disk in `ddt_dump/resuspension_mechanics/`, read 2026-08-29:

- **Warner, J. C., Sherwood, C. R., Signell, R. P., Harris, C. K. and Arango,
  H. G. (2008).** Development of a three-dimensional, regional, coupled wave,
  current, and sediment-transport model. *Computers & Geosciences* 34, 1284–1306.
  `warner_et_al_2008_sed_romscoastwst.pdf` — **vetted**
- **Sherwood, C. R., Aretxabaleta, A. L., Harris, C. K., Rinehimer, J. P.,
  Verney, R. and Ferré, B. (2018).** Cohesive and mixed sediment in the Regional
  Ocean Modeling System (ROMS v3.6) implemented in COAWST. *Geosci. Model Dev.*
  11, 1849–1871. — **vetted**

The questions they were compiled to answer, and the answers:

| Question | Answer |
|---|---|
| Which cases are idealized enough to adapt to particles? | Warner Ex. 1 (open-channel), Ex. 2 (trench, semi-idealized); Sherwood §3.2.1 (double resuspension), §3.2.2 (mixed bed), Fig. 2 (closure only) |
| Are they specified well enough to reproduce? | Yes — Warner Tables 2 and 4; Sherwood §3.2.1 gives every class property and the stress history |
| What can a Lagrangian model *not* reproduce? | Bed stratigraphy and active-layer thickness (no Eulerian bed); morphological feedback (no bed elevation update) |
| What is the right comparison instead? | Resuspension event count and intensity per class; deposition ordering; suspended-load partition between classes |

Standing rule this gate applied: **the parent study's observed sediment DDX
(Wu et al. 2025) is not a physics benchmark** and is excluded. It confounds
transport physics with source history and cannot isolate a physics defect.

Analytic references used for Tier A are in [REFERENCES.md](REFERENCES.md) and are
**unvetted** — Rouse, the well-mixed condition literature, Krone, van Rijn.

## 4. Design

### The two adaptation problems, stated up front

Comparing a particle model against Eulerian ROMS output is not free. Both of
these are known now rather than discovered mid-campaign:

1. **Particles need mass.** ROMS reports suspended-sediment concentration in
   kg m^-3. A cloud of unweighted numerical particles is not a concentration. Give
   each element a `mass` (M4 needs it anyway) and reduce with a mass-weighted
   kernel; or normalise every comparison and give up on absolute intensity. The
   first is better and cheap.
2. **ROMS erodes with a flux, we erode with a threshold.** ROMS uses
   `E = E_0 (1 - phi) (tau_b/tau_ce - 1)`, a *rate*. Our `resuspension()` lifts
   **every** settled element whose `tau_crit` is exceeded, every step — an
   all-or-nothing crossing. These do not produce the same event intensity, and
   the difference is not a tuning constant. See §7.

### Tiers

**Tier A — analytic.** No external data, seconds to run, CI.

| Case | Check |
|---|---|
| Rouse equilibrium profile | `C/C_a = [((h-z)/z)(a/(h-a))]^P` over a sweep of `P`; recover the exponent |
| Well-mixed condition | uniform stays uniform under strongly non-uniform `K(z)` with `w_s = 0` |
| Settling column | settled fraction equals `min(1, w_s t/h)` |
| BBL limits | SG2000 -> log law as `Ub -> 0`; -> pure wave as `Ur -> 0`; the three guard routes |
| Shields curve | `tau_crit_shields` against the Soulsby–Whitehouse fit |

The Rouse case is the single most diagnostic missing test: it exercises settling,
vertical mixing and resuspension height *jointly*, and `rouse.py` already assumes
exactly this profile, so a failure indicts the assumption itself.

**Tier B — published cases, adapted.**

*S3 — Sherwood Fig. 2 (closure check, no simulation; do first).* Cohesive
behaviour parameter `Pc(f_c)`; effective `tau_ce` per class for a bulk
`tau_cb = 0.1 Pa`; normalized excess-stress flux at `tau_b = 0.12 Pa`. Checks our
`mixed` `tau_crit` mode and its Jacobs/Yao ramp directly. Cheapest, highest value.

*S1 — Sherwood §3.2.1, the double resuspension experiment. Primary target.*
Fully specified: 20 m depth, 1DV (periodic lateral boundaries, flat bottom), no
floc dynamics; two stress events ~1.5 days apart lasting 1.5 and 1 days, the
first peaking at `tau_b = 1 Pa`; four classes —

| d (µm) | 4 | 30 | 62.5 | 140 |
|---|---|---|---|---|
| `tau_ce` (Pa) | 0.05 | 0.05 | 0.1 | 0.1 |
| `w_s` (mm/s) | 0.1 | 0.6 | 2 | 8 |

initial bed 41 x 1 mm layers, 25 % of each class. Metrics:

- **event number** — two events; the second resuspends only minimal 140 µm sand;
- **intensity** — fines dominate suspended load, only a small fraction of the
  coarsest class;
- **termination and settling** — on stress subsiding, coarse deposits first and
  fines stay up, giving fining-upward grading; the 4 µm class is still mostly
  suspended at day 5; net erosion 5 mm;
- **lag layer** — the two coarsest classes resist erosion (emergent armoring).

*Adaptation:* ROMS's stratigraphy and active-layer thickness become the
**deposition-order / burial-depth ordering of settled particles**. This is why S1
is simultaneously M1's sharpest test and the concrete argument for M2's
`burial_depth` state — the diagnostic the case demands is the state M2 adds.

*S2 — Sherwood §3.2.2, mixed-bed single stress event.* Exercises the cohesive and
`mixed` `tau_crit` paths.

*W1 — Warner Ex. 1, steady uniform open-channel flow.* Table 2: 10 000 x 100 x
10 m, `S_0 = 4e-5`, `u_bar = 1 m/s`, `z_0 = 0.0053 m`, `w_s = 1 mm/s`,
`tau_ce = 0.05 N/m^2`, `E_0 = 5e-5 kg/m^2/s`, porosity 0.90. Independently
`u* = sqrt(g h S_0) = 0.0626 m/s`, matching his Table 3 (0.0626 for k-eps) — so
the case self-checks. Note the implied Rouse number is small,
`P = w_s/(kappa u*) ~ 0.04`, i.e. a nearly uniform profile: this is a *weak*
discriminator of the exponent and a *strong* check that the model does not
spuriously stratify. Pair it with Tier A's `P` sweep, do not rely on it alone.

*W2 — Warner Ex. 2, migrating trench (van Rijn 1987 flume).* 30 m channel,
`D50 = 140 µm` mobile sand, `u_bar = 0.51 m/s`. **We cannot do morphology** — no
bed-elevation feedback — so the comparison is restricted to suspended-sediment
profiles at the five measurement stations and the erosion/deposition *pattern*
(accumulation at the upstream end, erosion at the downstream end). The limitation
is stated, not fudged.

**Tier C — keep green.** SG2000 Table 2 and settling vs Maggi (2013), already in
`tests/`.

### Harness

`tests/models/sppm_helpers.py` already fabricates models and environments without
a reader; M1 extends it with a 1DV column driver: constant and oscillating
currents, a parabolic / k-eps eddy diffusivity, a prescribed `tau_b(t)`, and
synthetic `Hs`/`Tp` so the SG2000 wave path is exercised without a wave reader.
This also retires the two unported tests that needed real MITgcm forcing.

## 5. Implementation steps

1. Extend the harness to a 1DV column driver.
2. S3 (closure only) — fastest signal on the cohesive/mixed threshold path.
3. Tier A, in order: Shields, settling column, well-mixed, Rouse sweep, BBL limits.
4. Per-element `mass` and a mass-weighted concentration reduction.
5. S1, the primary target; then S2.
6. W1, then W2.
7. `docs/sppmdrift/validation/REPORT.md` and the verdict.

## 6. Validation

Tier A and S3 run in CI (seconds). S1, S2, W1, W2 carry the `slow` marker.
Figures reproduce the corresponding published panel beside model output.
`REPORT.md` tabulates case, metric, tolerance, published value, model value,
pass/fail.

## 7. Open questions & decisions log

- **The pickup-flux question (likely to force a decision).** Matching *event
  intensity* against an Eulerian erosion rate may be impossible with an
  all-or-nothing threshold crossing. The Ariathurai–Partheniades stub now works
  (`calc_upward_resuspension_velocity`, M0), but its call site is still commented
  out and its normalisation `(tau - tau_c)/tau` differs from ROMS's
  `(tau/tau_c - 1)`. If S1 cannot match event intensity, pulling probabilistic
  erosion forward from M2/B4 into M1 is the likely resolution. **Not committed.**
- **The live blocker M1 must speak to.** `ddt_dump/RECAP.md` records that the two
  fastest settling classes lose ~45 % of their particles out of the domain to
  resuspension-driven hopping, which currently blocks the DDT class-weight
  inversion. S1 is the direct diagnostic: over-resuspending the coarse class
  relative to Sherwood's second event is the same defect. Report on it explicitly.
- **Which stress does the threshold see?** `bbl_stress` defaults to `combined`
  (the peak `rho u*cw^2`), but `tau_crit` was calibrated against a mean stress.
  `BBL_SG2000_PLAN.md` §9 lists this as still open. Tier B may settle it.
- **`cohesive_strength_coeff = 100` is uncalibrated** and sets the absolute
  cohesive threshold. S3 and S2 are the first real constraint on it. Calibration
  itself is M2.
- **RESOLVED by the campaign:** the pickup-flux question above is now answered —
  F2 and F3 in the report show a threshold cannot reproduce event intensity or a
  plausible lift-off count, and they share a cause. Probabilistic erosion is
  **required**, not optional. It remains M2/B4 work, but M2 should treat it as the
  module's first item rather than an optional extra.
- **NEW, from the campaign:** the reference model has a critical shear stress for
  *deposition* (`tau_d`, Krone 1962) that we lack entirely (F4).
- **Decided:** Wu et al. (2025) excluded from physics validation.
- **Decided:** synthetic wave forcing for M1; no WW3 reader.

## 8. Status

**Tier A and Tier B cases S3 and S1 built and run, 2026-08-29.** Full results and
findings in [validation/REPORT.md](validation/REPORT.md); one illustrated
notebook per case in [validation/notebooks/](validation/notebooks/). Headline:

*The settling / mixing / deposition chain is sound* — Rouse exponent recovered,
well-mixed condition satisfied, analytic deposition rate matched, and Sherwood's
double-resuspension experiment reproduced for suspended-load ordering,
deposition sequence and retention of the finest class.

*The erosion side is not.* Four defects isolated:

| | Finding | Fix belongs to |
|---|---|---|
| F1 | `mixed` threshold is not Sherwood Eq. 6: cohesive term uses the particle's own diameter rather than the bulk bed stress, and the `max(..., tau_c)` floor is missing. 0.6x to 14x wrong. | M2/B3 |
| F2 | An all-or-nothing threshold cannot reproduce event-intensity scaling; a weaker second event resuspends 0.92x as much sand as the first, where the reference gives "minimal". | M2/B4 |
| F3 | Settled elements above threshold re-lift **every step** — 299 lift-offs per sand particle over two events. **This is the mechanism behind the parent study's ~45 % fast-class domain loss.** | M2/B4 |
| F4 | No critical shear stress for deposition (`tau_d`, Krone); we deposit unconditionally. Not anticipated in this plan. | M2 (new) |
| F6 | `times_resuspended` was `uint8` and wrapped silently at 255, corrupting its own diagnostic — and it is exported by the DDT production runs. Fixed to `uint32`. | **fixed** |

Remaining to build: S2 (blocked in practice by F1), W1, W2, and an idealized
wave-forced case.
