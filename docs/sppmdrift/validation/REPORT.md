# M1 validation report

Status 2026-08-29. Tier A (analytic) and Tier B cases S3 and S1 are built and
run; S2, W1 and W2 are not yet built.

**Illustrated form:** [`notebooks/`](notebooks/) has one notebook per case, with
the setup, the forcing, the result plotted against the reference, and a verdict.
Start at [`notebooks/00_overview.ipynb`](notebooks/00_overview.ipynb).

Run the suite with:

```
/home/smullersoares/anaconda3/envs/opendrift_dev/bin/python -m pytest \
    tests/models/test_sppm_analytic.py tests/models/test_sppm_roms_cases.py \
    -o addopts="" -q --run-slow
```

## Verdict

**The settling / vertical-mixing / deposition chain is sound.** It reproduces the
Rouse equilibrium profile, satisfies the well-mixed condition, deposits at the
analytic rate, and — in Sherwood's double-resuspension experiment — reproduces
the class-ordered suspended load, the fining-upward deposition sequence and the
retention of the finest class after five days.

**The erosion side is not.** Three distinct defects were isolated, all of them
consequences of resuspension being an all-or-nothing threshold crossing rather
than a flux, plus one formulation error in the mixed-bed threshold. One of them
is the direct mechanism behind a blocker in the parent DDT study.

## Results

| Tier | Case | Metric | Tolerance | Reference | Model | Verdict |
|---|---|---|---|---|---|---|
| A | Parabolic diffusivity | shape, peak value | 1 % | `kappa u* h/4` | matches | **pass** |
| A | Settling column | settled fraction after 2 h | 0.05 | 0.72 (`w_s t/h`) | within tol | **pass** |
| A | Well-mixed condition | max bin excursion, two regions | 20 % | uniform | 7.5 % / 10 % | **pass** |
| A | Rouse profile, P = 0.25 | fitted exponent | 0.55–1.15 P | 0.25 | 0.223 | **pass** |
| A | Rouse profile, P = 0.50 | fitted exponent | 0.55–1.15 P | 0.50 | 0.404 | **pass** |
| B | S3 low-mud limit | tau_ce vs Shields | 1e-5 rel | identical | identical | **pass** |
| B | S3 monotonicity in mud | tau_ce non-decreasing | — | — | holds | **pass** |
| B | S3 cohesive limit vs Eq. 6 | tau_ce ratio | 25 % | max(tau_c, tau_cb) | 0.6–14 x | **FAIL (F1)** |
| B | S1 particle conservation | elements lost | 0 | — | 0 | **pass** |
| B | S1 suspended load ordering | peak fraction vs size | monotonic | fines dominate | 1.00/1.00/0.97/0.66 | **pass** |
| B | S1 deposition order | order of clearing | by `w_s` | coarse first | 42.0/45.5/61.0/never h | **pass** |
| B | S1 fines retained at day 5 | suspended fraction | > 0.5 | "mostly suspended" | 0.64 | **pass** |
| B | S1 second-event suppression | peak2 / peak1, 140 um | < 0.25 | "minimal" | 0.92 | **FAIL (F2)** |
| B | S1 lift-offs per particle | mean, 140 um | a few | a few | 299 | **FAIL (F3)** |
| C | SG2000 Table 2 | max abs error, 15 entries | 0.15 | Table 2 | 0.000 | **pass** |
| C | Settling vs Maggi (2013) | log10-RMSE, bb16 csf 0.7 | < 0.25 | — | 0.202 | **pass** |
| C | Settling vs Maggi (2013) | log10-RMSE, dietrich | < 0.25 | — | 0.219 | **pass** |

**A note on the well-mixed test.** As first written it checked only the interior
10-90 % of the column, which quietly excluded almost all the variation in `K` --
that band spans only a ~3x range. It now also checks the top 30 %, where the
surface is a genuine reflecting boundary and `K` falls by more than four orders of
magnitude; that is where a scheme missing the `dK/dz` correction would fail. The
model passes both (7.5 % and 10 % maximum bin excursion).

## Findings

### F1 — the `mixed` threshold is not Sherwood's mixed-bed formulation

*Notebook:* [`notebooks/04_S3_mixed_bed_closure.ipynb`](notebooks/04_S3_mixed_bed_closure.ipynb)

Sherwood et al. (2018) Eq. 6 is

```
tau_ce = max[ Pc * tau_cb + (1 - Pc) * tau_c , tau_c ]
```

where `tau_cb` is the **bulk critical stress of the bed** — one value for every
size class — and the result is floored at the particle's own `tau_c`.

Ours (`sedimentdrift.py:768`, mode `mixed`) is

```
tau_crit = (1 - Pc) * tau_shields(d) + Pc * tau_floc(d, Df, phi)
```

Both terms are per-particle, and `tau_floc` scales with the particle's **own**
diameter, so a sand grain embedded in a mud matrix is assigned the fictitious
strength of a 140 um floc rather than the strength of the mud holding it. There
is also no floor at `tau_c`. In the fully cohesive limit, against Eq. 6 with
`tau_cb = 0.1 Pa`:

| class | 4 um | 30 um | 62.5 um | 140 um |
|---|---|---|---|---|
| Sherwood Eq. 6 (Pa) | 0.100 | 0.100 | 0.120 | 0.159 |
| ours (Pa) | 0.064 | 0.478 | 0.995 | 2.229 |
| ratio | **0.6** | 4.8 | 8.3 | **14.1** |

So in a muddy bed we make fines *too easy* to erode (no floor) and coarse grains
**14x too hard**. Fix belongs to M2/B3.

**Scope.** `mixed` is not the production default (`auto` is), so the DDT
production runs are not affected through this path. But `auto` uses the same
`tau_crit_floc(d)` for every element it classes as cohesive, i.e. everything
below the 63 um cutoff — which is how a 31 um quartz grain in those runs acquires
a 0.48 Pa threshold against a legacy value of 0.09 Pa. That is consistent with
the near-zero resuspension already observed there.

### F2 — a threshold cannot reproduce event-intensity scaling

*Notebook:* [`notebooks/05_S1_double_resuspension.ipynb`](notebooks/05_S1_double_resuspension.ipynb)

Sherwood's second, weaker stress pulse "only resuspended minimal amounts of the
140 um sand", because erosion there is a flux,
`E = E_0 (1 - phi) (tau_b/tau_ce - 1)`, which scales with the excess stress.

Ours lifts every settled element of a class as soon as `tau_b > tau_ce`, however
small the excess. At the second event's 0.45 Pa peak against a 0.10 Pa threshold
we resuspend 0.92x as much sand as the 1.0 Pa first event — essentially no
suppression at all.

This is the adaptation problem anticipated in M1 §4, and it is not a tuning
constant: **matching event intensity requires a probabilistic pickup flux**
(M2/B4). The stub now works (`calc_upward_resuspension_velocity`) but is not
wired into `resuspension()`, and its excess-stress normalisation
`(tau - tau_c)/tau` differs from ROMS's `(tau/tau_c - 1)`.

### F3 — settled elements re-lift on every time step

*Notebook:* [`notebooks/05_S1_double_resuspension.ipynb`](notebooks/05_S1_double_resuspension.ipynb)

Mean lift-offs per particle over the whole two-event experiment:

| class | 4 um | 30 um | 62.5 um | 140 um |
|---|---|---|---|---|
| lift-offs per particle | 2.8 | 10.6 | 40.8 | **298.8** |

Two stress events should give a few lift-offs per particle, not three hundred. A settled
element whose threshold is exceeded is resuspended *again on the next step*, and
the coarse classes — which fall back to the bed quickly — cycle continuously.

**This is the mechanism behind the parent study's blocker.** `ddt_dump/RECAP.md`
records that the two fastest settling classes lose ~45 % of their particles out
of the domain, which currently blocks the class-weight inversion. Each hop
displaces a particle downstream; three hundred hops is a large random-walk
displacement that has nothing to do with the physics. Isolated here in a case
with a known answer, the defect is unambiguous.

F2 and F3 share a cause and would share a fix.

### F6 — `times_resuspended` was an 8-bit counter and wrapped silently

Found while writing the notebooks, when the same quantity came out as 44 in one
analysis and 254 in another. `SedimentElement.times_resuspended` was declared
`np.uint8`, so it **wraps at 255**: the live element array showed a mean of 44.6
for the coarse class while the exported history peaked at exactly 255.0. Neither
number was real — the true mean is **298.8**, with a maximum of 341.

The counter is now `np.uint32` (`sedimentdrift.py:146`). Everything reported for
F3 above uses the corrected values.

**This reaches beyond the validation suite.** `times_resuspended` is in the
`export_variables` list of the DDT production runs, so any analysis of
resuspension counts in existing output is wrapped wherever a particle exceeded
255 lift-offs — which, given F3, the fast classes certainly did. Those outputs
cannot be re-derived without re-running, but the affected diagnostic should not
be trusted.

### F4 — no critical shear stress for deposition

Sherwood implements an optional `tau_d` (Krone 1962): deposition is zero while
`tau_b > tau_d` and increases linearly as `tau_b` falls. We deposit
unconditionally whenever an element reaches the bed. The effect of `tau_d` is to
keep material in suspension in the bottom layer during an event. This was not
anticipated in the M1 plan and is a new item for M2.

### F5 — the Rouse recovery is biased by the settle-and-resuspend cycle

*Notebook:* [`notebooks/03_rouse_profile.ipynb`](notebooks/03_rouse_profile.ipynb)

The recovered exponent is biased low (0.223 for P = 0.25; 0.404 for P = 0.50)
because elements resting on the bed awaiting the next resuspension are not part
of the suspended profile. A sweep over the treatment of that population:

| P | outer dt | duration | bed-resting elements | fitted P |
|---|---|---|---|---|
| 0.25 | 300 s | 24 h | 267 | 0.223 (excluded) |
| 0.25 | 300 s | 24 h | 267 | 0.582 (counted at z_a) |
| 0.25 | 60 s | 24 h | 54 | 0.407 (counted at z_a) |
| 0.25 | 60 s | 48 h | 47 | 0.443 (counted at z_a) |

Counting them at the reference height over-corrects badly; excluding them is much
the better estimator, and is what the shipped test does. Shrinking the time step
shrinks the bed-resting population as expected (267 → 54 for a 5x smaller step),
confirming it is a discretisation artifact of the lift-off cycle rather than a
physical result. At P = 1 the suspended population becomes too small for a
meaningful fit, so the benchmark is scoped to P <~ 0.5.

## Not yet built

- **S2** — Sherwood section 3.2.2, mixed-bed single stress event. Blocked in
  practice by F1: comparing a mixed-bed case is not meaningful until the
  threshold formulation is fixed.
- **W1** — Warner Example 1, steady open-channel flow. Tier A's Rouse sweep
  already covers the same balance; W1 adds a published profile and a
  turbulence-closure sensitivity.
- **W2** — Warner Example 2, migrating trench. Restricted to suspended-sediment
  profiles and the erosion/deposition pattern; morphology is out of reach.
- **Wave path** — the SG2000 wave branch is exercised only by unit tests
  (`test_sedimentdrift_bbl.py`); no idealized wave-forced case yet.
