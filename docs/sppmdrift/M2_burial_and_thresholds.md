# M2 — Burial and resuspension thresholds

## 1. Goal & scope

Two halves, both of which make settled material behave like a *deposit* rather
than like a particle parked on a surface:

1. **Burial** — a simple scheme that makes material settled for a long time
   harder, and eventually impossible, to resuspend.
2. **Thresholds** — improve the physics and the *environmental* setting of the
   resuspension threshold, which today is driven by scalar config defaults.

Design principle, inherited: **stay Lagrangian; no Eulerian bed grid.** The
existing work deliberately avoided a layered bed as "the expensive option we will
avoid", and that choice should be preserved unless the literature gate overturns
it. "As simple as possible" is a requirement, not an aspiration.

Out of scope: a full stratigraphic bed model; diagenesis.

## 2. Current state

**Burial does not exist.** `burial` appears once in the whole DDT project, and
rhetorically — `RESUSPENSION_MECHANICS_PLAN.md` L110, "the 'fresh DDT mud erodes
easier than buried mud' effect". `consolidat` and `armor` have zero hits in the
project's own code.

The one adjacent hook is consolidation:
`floc_strength.consolidate_phi(phi0, t_settled, t_consol, phi_max)` relaxes the
solids fraction from `phi0` toward `phi_max`, and `_critical_shear_stress` feeds
that into the cohesive fractal strength law, so a deposit does stiffen with age.
But **nothing is ever removed from the resuspendable pool**: an element settled
for a decade is still a candidate for lift-off every step.

Thresholds today (`sedimentdrift.py:768`), all reading scalar config:

| Config | Default | Comment |
|---|---|---|
| `tau_crit_mode` | `auto` | routes per element by `sed_class` / 63 µm cutoff |
| `cohesive_strength_coeff` | 100.0 | **admittedly O(100)-arbitrary; sets the absolute cohesive threshold** |
| `consolidation_timescale` | 86400 s | |
| `consolidated_solids_fraction` | 0.4 | |
| `bed_mud_fraction` | 0.5 | a single number for the whole domain |
| `mud_fraction_lower` / `upper` | 0.05 / 0.35 | Jacobs (2011) / Yao (2022) ramp |

`vertical_mixing:erosion_rate_constant` and `bed_porosity` now exist (M0) and
feed the currently-unused Ariathurai–Partheniades pickup stub.

## 3. Literature gate — OPEN

**Blocking. Compile and review jointly before any design.**

Questions this gate must answer:

1. What is the simplest defensible representation of burial for a Lagrangian
   model — a hard mixed-layer cut, or a smooth exposure probability? What does
   the layered-bed literature actually require that we can legitimately drop?
2. Should burial be *reversible*? Over the 70-year timescales this model is aimed
   at, bioturbation demonstrably re-exposes buried material. What sets the rate?
3. Should consolidation be driven by overburden (`burial_depth`) rather than by
   elapsed time (`time_since_settled`)? The physical driver is effective stress.
4. What actually calibrates `cohesive_strength_coeff`? Which flume datasets are
   usable given the floc properties we seed?
5. Is armoring something to parameterise, or should it emerge from a per-element
   grain-size population plus burial? (Preference: emergent, then *check* it.)
6. Threshold or rate? Does a per-particle threshold crossing need replacing by a
   probabilistic pickup flux — and does M1 already force that?

Candidate papers — **all unvetted**; bibliographic details to be confirmed at the
gate. See [REFERENCES.md](REFERENCES.md).

*Layered beds, consolidation, armoring*
- Sanford (2008), modelling a dynamically varying mixed sediment bed with erosion,
  deposition, bioturbation, consolidation and armoring — **the central reference**
- Sanford & Maa (2001), a unified erosion formulation for fine sediments *(HTML copy on disk)*
- Harris & Wiberg (2001), suspended transport and bed reworking on shelves
- Warner et al. (2008) §bed model; Sherwood et al. (2018) revised stratigraphy *(both on disk)*
- Le Hir et al. (2011), erosion of cohesive/mixed sediment
- Winterwerp & van Kesteren (2004), physics of cohesive sediment *(book)*

*Bioturbation and mixed-layer depth*
- Boudreau (1998), mean mixed depth of sediments — the ~9.8 cm global mean
- Boudreau (1986), diffusive models of bioturbation
- Meysman et al. (2006); Teal et al. (2008), global bioturbation intensity
- Wheatcroft et al. (1990), tracer burial

*Thresholds, hiding/exposure, mixed beds*
- Wiberg & Smith (1987) *(on disk)*; Soulsby & Whitehouse (1997)
- Jacobs et al. (2011); Yao et al. (2022) *(both on disk)*
- Winterwerp et al. (2012) *(on disk)*; Kranenburg (1994)
- Wiberg, Drake & Cacchione (1994), resuspension and bed armoring under high stress

*Erosion rate*
- Partheniades (1965); Ariathurai & Arulanandan (1978) *(on disk)*; Krone (1962)

*Pollutant-specific (why reversibility matters here)*
- Eganhouse & Pontolillo (2008), DDE in Palos Verdes shelf sediments

## 4. Design

*Stub — pending §3.* Sketch only, to be confirmed, refuted or replaced:

- **B0** `burial_depth` per element, growing while settled at the local net
  deposition rate (config or 2D map first; self-consistent from the model's own
  deposition flux later).
- **B1** burial gates resuspension: a soft exposure probability decaying with
  `burial_depth / L_mix`, with a hard cut as a `simple` mode. This is the
  Lagrangian analogue of Sanford's `tau_c(z)` bed profile without the bed grid.
- **B2** bioturbational re-exposure as a random walk on `burial_depth`.
- **B3** environmental thresholds: `bed_mud_fraction` promoted from scalar config
  to a reader variable; `cohesive_strength_coeff` calibrated; hiding/exposure.
- **B4** optional probabilistic erosion (Ariathurai–Partheniades pickup flux).

## 5. Implementation steps

*Stub — pending §3.*

## 6. Validation

Sherwood §3.2.1 (M1's S1) is *already* the natural test: its fining-upward
grading, lag layer and second-event suppression are exactly burial and armoring
signatures, and M1 can only compare them by proxy without `burial_depth`. Sherwood
§3.2.3 adds biodiffusion variants. Expect M1 to hand M2 a concrete failure to fix.

## 7. Open questions & decisions log

- Does a soft exposure probability introduce a spurious timescale that a hard cut
  avoids? Both are approximations to a `tau_c(z)` profile.
- Deposition rate from a map (defensible, measured) vs from the model's own flux
  (self-consistent, but couples particles and needs M4's mass).
- Whether B4 gets pulled forward into M1 (see M1 §7).

## 8. Status

**Planned. Literature gate open — do not start §4.**
