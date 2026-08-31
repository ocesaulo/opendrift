# M2 — Bed physics: erosion flux, thresholds, and burial

> **Revised 2026-08-31, after the M1 validation campaign.** M1 found that the
> erosion step is not merely miscalibrated but **non-convergent** — the amount
> resuspended is set by the timestep, not by the flow. That changes this module's
> shape: what was an optional last item (a pickup flux) is now the prerequisite
> for everything else in it, and one new item (a deposition threshold) was added.
> §7 records a **scope decision for review** before work starts.

## 1. Goal & scope

Make settled material behave like a **deposit** rather than like a particle
parked on a surface. Three parts, in the order they now have to be done:

1. **Erosion as a flux.** Replace the per-particle threshold crossing with a
   rate, so that how much material leaves the bed depends on the flow and not on
   the timestep.
2. **Thresholds.** Fix the mixed-bed formulation, give the threshold a real
   environmental setting, and calibrate the one arbitrary constant in it.
3. **Burial.** A simple scheme making material settled for a long time harder,
   and eventually impossible, to resuspend — the module's original goal.

**Why that order.** Burial gates resuspension; if resuspension itself does not
converge, gating it is gating a numerical artifact. The same applies to threshold
work: a better threshold changes *when* the crossing happens, not the fact that
crossing it mobilises an arbitrary amount of material. Part 1 is therefore
load-bearing for parts 2 and 3, not merely more urgent than them.

Design principles, inherited and still binding:

- **Stay Lagrangian.** No stratigraphic bed, no layered `tau_c(z)`, no bed-cell
  state machine. The existing work deliberately avoided a layered bed as "the
  expensive option we will avoid". Part 1 needs *one* Eulerian-ish quantity — an
  erodible mass per unit area — and that is the only concession; see §4.
- **As simple as possible**, and that is a requirement rather than an aspiration.
- **`SedimentDrift` stays behaviour-identical by default.** Every change lands
  behind a config whose default reproduces today's behaviour, so the DDT
  production brackets remain reproducible.

Out of scope: a full stratigraphic bed model; diagenesis; morphological feedback
(no bed-elevation update).

## 2. Current state

### What M1 established

M1 is the reason this document changed. Full results in
[validation/REPORT.md](validation/REPORT.md); narrative assessment in
[validation/first_claude_assessment_M1_physics.txt](validation/first_claude_assessment_M1_physics.txt).

| | Finding | Lands in |
|---|---|---|
| **F3b** | **The lift-off rate does not converge under timestep refinement**: +63 % over a 4× refinement of `dt`, monotonic, no limit. The count is set by how often `resuspension()` is called. | **B1** |
| F3 | Settled elements above threshold re-lift *every step* — 299 lift-offs per coarse particle over two events. The mechanism behind the parent study's ~45 % fast-class domain loss, which is therefore a numerical artifact. | **B1** |
| F2 | No intensity scaling: a 0.45 Pa event resuspends 0.92× as much sand as a 1.00 Pa event, where the reference gives "minimal amounts". | **B1** |
| F1 | The `mixed` threshold is not Sherwood Eq. 6 — the cohesive term uses the particle's own diameter instead of the bulk bed stress, and the `max(·, τ_c)` floor is missing. 0.6× to 14× wrong. | **B2** |
| F4 | No critical shear stress for *deposition* (Krone `τ_d`); we deposit unconditionally. | **B3** |
| F6 | `times_resuspended` was `uint8` and wrapped at 255. | **fixed** |

What M1 also established, and which this module can rely on: settling, vertical
mixing, the settling⊗mixing balance and the SG2000 bed stress are all validated
against external references. **The inputs to erosion are trustworthy; the erosion
step is not.**

### The code as it stands

`resuspension()` (`sedimentdrift.py:841`) advances the consolidation clock, gets a
threshold from `_critical_shear_stress()` and a stress from `calc_bottom_stress()`,
and lifts **every** element where `bottom_stress > threshold` — unconditionally,
once per outer timestep. There is no rate, no erodible-mass limit, and no
dependence on the size of the excess.

**Burial does not exist.** `burial` appears once in the whole DDT project, and
rhetorically (`RESUSPENSION_MECHANICS_PLAN.md` L110, "the 'fresh DDT mud erodes
easier than buried mud' effect"). `consolidat` and `armor` have zero hits in the
project's own code.

The one adjacent hook is consolidation:
`floc_strength.consolidate_phi(phi0, t_settled, t_consol, phi_max)` relaxes the
solids fraction toward `phi_max`, and `_critical_shear_stress` feeds it into the
cohesive fractal strength law, so a deposit does stiffen with age. But **nothing
is ever removed from the resuspendable pool**: an element settled for a decade is
still a lift-off candidate on every step.

Thresholds today (`sedimentdrift.py:771`), all reading scalar config:

| Config | Default | Comment |
|---|---|---|
| `tau_crit_mode` | `auto` | routes per element by `sed_class` / 63 µm cutoff |
| `cohesive_strength_coeff` | 100.0 | **admittedly O(100)-arbitrary; sets the absolute cohesive threshold** |
| `consolidation_timescale` | 86400 s | |
| `consolidated_solids_fraction` | 0.4 | |
| `bed_mud_fraction` | 0.5 | a single number for the whole domain |
| `mud_fraction_lower` / `upper` | 0.05 / 0.35 | Jacobs (2011) / Yao (2022) ramp |

Available and unused: `vertical_mixing:erosion_rate_constant` (E₀) and
`vertical_mixing:bed_porosity`, added in M0, feeding
`calc_upward_resuspension_velocity` — a working but uncalled
Ariathurai–Partheniades stub whose excess-stress normalisation
`(τ − τ_c)/τ` differs from the ROMS/COAWST form `(τ/τ_c − 1)`.

**Elements carry no mass.** This is the non-obvious prerequisite for B1; see §4.

## 3. Literature gate — OPEN

**Blocking. Compile and review jointly before any design is fixed.** M1 answered
some of the questions this gate was originally posed to answer, so they are split
below into what is now settled by measurement and what remains open.

### Settled by M1 — do not relitigate

- **Threshold or rate?** Rate. F2, F3 and F3b are three faces of one defect, and
  F3b in particular is not a calibration problem: a threshold crossing has no
  convergent limit. Probabilistic erosion is **required, not optional**.
- **Is armoring worth parameterising?** No — treat it as an emergent diagnostic.
  M1's S1 already shows a coarse lag forming from differential settling alone.

### Still open — what the papers must answer

1. **Which excess-stress form?** ROMS/COAWST uses `E = E₀(1−φ)(τ_b/τ_c − 1)`; our
   stub uses `(τ − τ_c)/τ`. They differ substantially near threshold. Which is
   better supported, and does `E₀` carry the same meaning in each?
2. **What sets `E₀` for our material?** It is currently a config default of
   5e-5 kg m⁻² s⁻¹ with no provenance. Which flume datasets are usable given the
   floc properties we seed?
3. **Erodible mass and the active layer.** ROMS limits erosion to an active-layer
   thickness. What is the minimum defensible Lagrangian equivalent (see §4), and
   is an area-based accounting acceptable without a layered bed?
4. **Burial: hard cut or smooth exposure?** What does the layered-bed literature
   actually require that we can legitimately drop?
5. **Is burial reversible?** Over the 70-year timescales this model targets,
   bioturbation demonstrably re-exposes buried material. What sets the rate, and
   what mixed-layer depth?
6. **Should consolidation drive the *bulk bed* stress rather than each particle's
   own strength?** M1's F1 suggests the current arrangement conflates the two
   (see §4/B2). The physical driver is effective stress, i.e. overburden.
7. **What calibrates `cohesive_strength_coeff`?**
8. **Deposition threshold.** Is Krone's `τ_d` the right form, and is the
   "typically about half `τ_c`, but unrelated to it" rule of thumb defensible?

### Candidate papers — all `unvetted` unless noted

*Erosion rate / pickup flux (now the priority)*
- Ariathurai & Arulanandan (1978), erosion rates of cohesive soils *(on disk)*
- Partheniades (1965), erosion and deposition of cohesive soils
- Sanford & Maa (2001), a unified erosion formulation for fine sediments *(on disk, HTML)*
- Winterwerp et al. (2012) *(on disk)*; Mehta (2013), *Hydraulics of Fine Sediment Transport*
- Warner et al. (2008) §bed model; Sherwood et al. (2018) §2.3.2 *(both on disk, both read)*

*Deposition threshold*
- Krone (1962), flume studies of estuarial shoaling
- Whitehouse et al. (2000); Spearman & Manning (2008)

*Layered beds, consolidation, armoring*
- **Sanford (2008)** — the central reference for a dynamically varying mixed bed
- Harris & Wiberg (2001); Le Hir, Cayocca & Waeles (2011)
- Winterwerp & van Kesteren (2004), *Physics of Cohesive Sediment* (book)
- Wiberg, Drake & Cacchione (1994), armoring under high bottom stress

*Bioturbation and mixed-layer depth*
- Boudreau (1998), mean mixed depth of sediments — the ~9.8 cm global mean
- Boudreau (1986); Meysman, Middelburg & Heip (2006); Teal et al. (2008)

*Thresholds, hiding/exposure, mixed beds*
- Wiberg & Smith (1987) *(on disk)*; Soulsby & Whitehouse (1997)
- Jacobs et al. (2011); Yao et al. (2022) *(both on disk)*; Kranenburg (1994)

*Pollutant-specific (why reversibility matters here)*
- Eganhouse & Pontolillo (2008), DDE in Palos Verdes shelf sediments

## 4. Design

*Sketch, pending §3. The B1 formulation is worked through further than the rest
because M1 constrains it tightly; it is still a proposal.*

### B1 — erosion as a convergent flux  *(the priority)*

The pickup flux gives a mass per unit bed area per unit time:

```
E = E₀ (1 − φ) (τ_b/τ_c − 1)        [kg m⁻² s⁻¹]
```

Turning that into Lagrangian lift-off decisions:

```
M_available  = Σ mᵢ over settled elements in a bed patch of area A
M_eroded     = E · A · Δt
p            = min(1, M_eroded / M_available)     # per-element lift probability
lift if U(0,1) < p
```

**Why this converges.** `p ∝ Δt`, so the expected mass eroded per step is
`E·A·Δt` and over an event of duration `T` it is `E·A·T` — independent of `Δt`.
That is precisely the property F3b shows the current scheme lacks.

**Caveat, and it must be tested rather than assumed.** The argument holds only
while `p < 1`. Where the flow is energetic enough to saturate it (`τ_b ≫ τ_c`, or
too few particles carrying too little mass in a patch), every element lifts every
step and the timestep dependence returns — the scheme degrades exactly back to
today's behaviour. That may be defensible as "full mobilisation", but it is a
genuine limitation of the formulation, not a free pass. The convergence gate in
§6 should therefore sweep a stress range that spans both regimes and report where
saturation begins, rather than testing one comfortable operating point.

Two prerequisites, both real work:

- **Per-element `mass`.** Elements carry none today. M1 already needed it to
  compare a particle cloud against an Eulerian SSC, and **M4 needs it regardless**
  — so this pulls a piece of M4's state into M2. Flagged as a plan-level
  dependency change.
- **A bed area `A`.** Three options, with a recommendation:
  1. *A coarse binning grid used only for the erosion denominator.* Recommended.
     It is not a stratigraphic bed — no layers, no per-cell state carried between
     steps, just a denominator evaluated on the fly — so it does not breach the
     "stay Lagrangian" rule in spirit. Cost is one `histogram2d`-style reduction
     per step over settled elements, which `ddt_sim/analysis.py` already does
     vectorized for the deposition maps.
  2. *A fixed footprint per element.* Simplest, but makes the erosion rate depend
     on particle count, which is a numerical parameter — trading one artifact for
     another.
  3. *A kernel density estimate over neighbours.* Most defensible, most
     expensive, and shares machinery with M4/M5's colocation kernel — worth
     noting as the eventual convergence point rather than the first step.

Config: `vertical_mixing:erosion_mode = threshold (default, today's behaviour) |
flux`. Reuses the existing `erosion_rate_constant` and `bed_porosity`, and revives
`calc_upward_resuspension_velocity` after settling question 1 in §3.

### B2 — the threshold, corrected  *(F1)*

Adopt Sherwood Eq. 6's structure:

```
τ_ce = max[ P_c·τ_cb + (1 − P_c)·τ_c , τ_c ]
```

Two changes: the cohesive term becomes a **bulk bed** critical stress `τ_cb`,
identical for every size class, and the `max(·, τ_c)` floor is restored.

This suggests a unification worth testing rather than assuming: **consolidation
should raise `τ_cb`, not each particle's own floc strength.** That is more
physical (the driver is effective stress in the deposit, not the strength of an
individual floc) and it removes the conflation that produced F1. Whether `τ_cb`
comes from `consolidate_phi` + Kranenburg, from a config, or from a reader field
is a gate question.

Also here: promote `bed_mud_fraction` from a scalar config to a **reader
variable**, so the `mixed` mode can use a real sediment map (`data/sediment/` in
the parent project); calibrate `cohesive_strength_coeff`; add Wiberg–Smith
hiding/exposure for mixed grain sizes.

### B3 — deposition threshold  *(F4, new, small)*

Suppress deposition while `τ_b > τ_d`, linearly ramping as `τ_b → 0`
(Whitehouse's "linear depositional flux"), with a constant-flux option. Keeps
material in suspension near the bed during an event. Config
`vertical_mixing:tau_crit_deposition`, default disabled for back-compatibility.

### B4 — burial  *(the module's original goal)*

- `burial_depth` per element, growing while settled at the local net deposition
  rate — from a config or 2D map first (defensible here: Wu et al. 2025 measured
  accumulation rates), later diagnosed from the model's own deposition flux.
- Burial gates resuspension through a **soft exposure probability** decaying with
  `burial_depth / L_mix`, with a hard cut as the `simple` mode. This is the
  Lagrangian analogue of Sanford's `τ_c(z)` bed profile without the bed grid.
  `L_mix` default from Boudreau (1998), ≈9.8 cm.
- **Bioturbational re-exposure** as a random walk on `burial_depth`, making burial
  reversible — essential for a 70-year DDT problem where cores show buried
  material resurfacing.

Note that B4 composes cleanly with B1: an exposure probability multiplies the
lift-off probability, rather than gating a binary decision.

## 5. Implementation steps

1. Clear the §3 literature gate (jointly). Settle questions 1–3 at minimum before
   touching B1.
2. **B1a** — per-element `mass`, seeded and exported; mass-weighted concentration
   diagnostic. Verify M1's S1 numbers are unchanged (mass is inert until B1b).
3. **B1b** — erodible-mass accounting and the `flux` erosion mode, behind a config
   defaulting to today's behaviour.
4. **B1 gate** — invert `test_s1_lift_off_rate_does_not_converge` (§6). *Nothing
   else in this module proceeds until this passes.*
5. **B2** — corrected mixed-bed threshold; un-xfail
   `test_s3_matches_sherwood_eq6_in_cohesive_limit` and retire
   `test_s3_divergence_is_characterised`.
6. **B3** — deposition threshold.
7. **B4** — burial depth, exposure gating, bioturbation.
8. Re-run the full M1 suite and rebuild the validation notebooks; update
   `REPORT.md` and the scoreboard in `00_overview.ipynb`.

## 6. Validation — acceptance gates

These are concrete and already partly written. M1 left the tests in place
asserting the *defects*, precisely so that M2 has falsifiable targets.

| Gate | Test | Today | Required after |
|---|---|---|---|
| **Convergence** | `test_s1_lift_off_rate_does_not_converge` | asserts +63 % over a 4× `dt` refinement | **inverted**: eroded mass agrees across `dt` to within sampling noise, over a stress sweep spanning both the unsaturated and saturated regimes |
| Intensity scaling | `test_s1_second_event_does_not_suppress_coarse_resuspension` | asserts peak2/peak1 > 0.75 | inverted to the docstring's `< 0.25` |
| Hopping | `test_s1_coarse_classes_hop_repeatedly` | asserts > 100 lift-offs/particle | a handful per event |
| Mixed-bed threshold | `test_s3_matches_sherwood_eq6_in_cohesive_limit` | strict `xfail` | passes within 25 % |
| No regression | the whole M1 suite | 46 pass, 1 xfail | still passes; Tier A must be untouched |
| Back-compat | RECAP C4 sanity config | — | reproduces its baseline with `erosion_mode='threshold'` |

New cases this module should add:

- **S2** — Sherwood §3.2.2, the mixed-bed single stress event. Deliberately not
  built in M1 because comparing a mixed-bed case is meaningless until F1 is
  fixed; it becomes the natural B2 acceptance case.
- **Erodible-mass conservation** — the mass eroded over an event must not exceed
  the mass available, and must be independent of particle count.
- **A burial case** — Sherwood §3.2.3 adds biodiffusion variants.

## 7. Open questions & decisions log

- **SCOPE DECISION FOR REVIEW.** M1 added a prerequisite (B1, the erosion flux)
  that was not in the original five-module conception, where module 2 was burial
  plus threshold physics. Two options:
  - **(a) M2 absorbs it**, as written here. Keeps the bed physics in one place;
    makes M2 substantially larger than originally scoped, and defers burial —
    the module's stated goal — behind two other items.
  - **(b) Split it out as its own module** (say M1.5 or M2a, "erosion as a
    flux"), leaving M2 as burial + thresholds per the original plan.

  *Recommendation: (b).* B1 is a self-contained, independently testable change
  with a single unambiguous acceptance gate, and it is a prerequisite for two
  other modules' validation rather than only this one's. Splitting keeps each
  module's "done" crisp. **Deferred to S. Soares.**
- **Resolved by M1:** probabilistic erosion is required, not optional; armoring
  stays an emergent diagnostic.
- **New dependency:** B1 pulls per-element `mass` from M4 into this module. M4's
  phase-mass vector then extends it rather than introducing it.
- Does a soft exposure probability introduce a spurious timescale that a hard cut
  avoids? Both approximate a `τ_c(z)` profile.
- Deposition rate for `burial_depth` from a map (defensible, measured) versus
  from the model's own flux (self-consistent, but couples particles and needs the
  same mass B1 introduces).
- Once erosion is a flux, is the Rouse lift-off height still the right release
  distribution, or should release height couple to the pickup rate?

## 8. Status

**Planned; literature gate open — do not start §4.** Revised 2026-08-31 to
incorporate M1's findings. A scope decision (§7) is awaiting review.
