# M5 — Dynamic particle size

## 1. Goal & scope

Let particle size evolve during a run instead of being fixed at seeding, either
from statistics and environment alone, or from a parameterised treatment of
particle-to-particle interaction.

Flagged from the outset as **low priority and possibly infeasible**. This document
therefore opens with a go/no-go rather than a design.

## 2. Current state

Nothing, and dynamic flocculation was explicitly excluded by the existing work:
*"No floc dynamics: turbulent shear breaks/re-aggregates flocs (evolving size &
density) — wholly outside a fixed-grain ODE."*

But the **state is already there**. Every element carries `grain_diameter`, `d0`,
`fractal_dim` and `phi0`; `settling.py` already maps them to an effective density
and a settling velocity; `floc_strength.py` already maps them to an erosion
threshold. M5 does not need new state — it needs those variables to stop being
constant. That is a much smaller change than it first appears, and it is why the
module is worth writing down even at low priority.

## 3. Literature gate — OPEN

**Blocking.** Questions:

1. Is a single-class kinetic equation for a characteristic floc size enough for
   our purposes, or is a two-class / full population balance required?
2. What drives it — turbulent shear `G` and concentration `C`. Can we obtain `G`
   honestly from the model's own diffusivity or dissipation?
3. Is a fixed fractal dimension tenable once size evolves? (Links to M3 Q2.)
4. **Feasibility of branch B**: is a neighbour search over ~1e5 particles per step
   affordable, and is the Lagrangian particle density even a valid estimator of
   physical SPM concentration?

Candidate papers — **all unvetted**:

*Single-class kinetics*
- Winterwerp (1998, 2002), turbulence-induced flocculation of cohesive sediment
- Khelifa & Hill (2006), variable fractal dimension
- Strom & Keyvani (2011); Keyvani & Strom (2014), growth and breakup under cycled shear

*Population balance*
- Smoluchowski (1917); Jackson (1990); Burd & Jackson (2009) *(review)*
- Lee et al. (2011), two-class population balance
- Verney et al. (2011), floc population over a tidal cycle *(FLOCMOD; the case
  Sherwood et al. 2018 §3.1.1 validates against)*
- Maerz et al. (2011)

*Kernel machinery (shared with M4)*
- Benson & Bolster (2016); Schmidt et al. (2017)

## 4. Design

*Stub — pending §3.* Two branches, deliberately separated:

- **Branch A — statistics only.** A per-element kinetic ODE for floc size driven
  by turbulent shear and an *environment* SPM concentration. No particle–particle
  interaction, so it is genuinely achievable, and it immediately feeds
  `settling_dynamic` (M3) and `fractal_dim`/`phi0` (M2).
- **Branch B — particle–particle.** Local concentration and collision statistics
  estimated from the Lagrangian cloud itself via a colocation kernel, then
  aggregation and breakup between nearby pairs.

**Structural dependency, recorded up front: branch B strictly depends on M4.** A
fixed number of numerical particles per release is *not* an estimator of physical
SPM concentration. It becomes one only when elements carry mass and the kernel is
mass-weighted — which is exactly what M4 provides. Attempting branch B before M4
would produce a concentration field with no physical units and no meaning.

Branch A has no such dependency and can be done any time after M3.

## 5. Implementation steps

*Stub — pending §3 and the go/no-go.*

## 6. Validation

Sherwood et al. (2018) §3.1.1 is a ready-made target and is already on disk: it
reproduces the Verney et al. (2011) laboratory tidal-shear experiment as a 0-D
case — 15 cohesive classes log-spaced 4–1500 µm, `nf = 1.9`, settling switched
off, constant SSC 0.093 kg m^-3, all mass initially in the 120 µm class, `G(t)`
cycling 0–12 s^-1 — and reports a 24 µm RMS difference in mass-weighted mean
diameter. A second case from Keyvani & Strom (2014) cycles `G = 15 s^-1` growth
against `G = 400 s^-1` disaggregation.

A 0-D case with settling off is about the friendliest possible test for a
Lagrangian model: it removes transport entirely and tests only the kinetics.
**This is the go/no-go experiment** — if branch A cannot reproduce it, the module
stops there.

## 7. Open questions & decisions log

- Cost of a per-step neighbour search at production particle counts (branch B).
- Whether `G` from a `constant`-diffusivity run is meaningful at all; the DDT runs
  use `diffusivitymodel='constant'`, which would give a constant floc size.
- Whether floc breakup on resuspension should be represented — an element lifted
  at 1 Pa plausibly does not survive intact.

## 8. Status

**Proposed. Go/no-go pending; literature gate open.** Branch A only after M3;
branch B only after M4.
