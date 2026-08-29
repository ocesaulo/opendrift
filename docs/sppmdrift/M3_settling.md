# M3 — Settling velocity

## 1. Goal & scope

Settling rates for the full range of particle types this model must carry —
cohesive aggregates that flocculate, and solid grains and debris that do not —
parameterised by **density, size, cohesiveness (a flag) and shape (for solids)**,
usable **offline as run setup** and eventually **dynamically** during a run.

Out of scope: the *kinetics* of aggregation, which is M5. M3 supplies the
velocity given a particle's state; M5 evolves that state.

## 2. Current state

The most advanced of the five modules. `opendrift/models/settling.py` provides
drag laws (Schiller–Naumann, Dietrich, BB16), an iterative vectorized
`terminal_velocity`, BB16 shape factors from a Corey shape factor, Maggi fractal
effective density with a permeable (Brinkman) variant, and a
`settling_velocity(model, ...)` dispatcher over
`{stokes, dietrich, bb16, maggi, maggi_permeable}`. Wired in through
`SedimentDrift.update_terminal_velocity` (`:685`) with configs
`vertical_mixing:settling_model` (default `prescribed`) and `settling_dynamic`
(default `False`). Validated against Maggi (2013) measurements at log10-RMSE 0.20,
with BB16's shape correction measurably improving on a sphere.

Known gaps, from `SETTLING_IMPLEMENTATION.md` §8 and the M0 audit:

- **Per-element model dispatch is deferred.** A run picks *one* closure. For a
  mixed population — flocs plus solid barrel debris and paint chips — that is a
  correctness limit, not a convenience one.
- **`settling_dynamic=True` has never been enabled or validated.** M5 needs it.
- **Cohesiveness is inferred, not declared**: from `sed_class`, or from a 63 µm
  size cutoff. The same property should route both the settling closure and the
  erosion threshold; today they are decided separately.
- **Shape is only exposed via the Corey factor.** `bb16_drag_factors` accepts
  full L/I/S axes but nothing plumbs them through.
- **Missing closures**: no Ferguson–Church, no Winterwerp floc settling, no
  Khelifa–Hill variable fractal dimension. `fractal_dim` is a fixed input.
- **The offline reference code is buggy.** `ddt_dump/codes/maggi.py` has a
  `D0_from_Df` that returns its own argument and a `maggi_H_approx` that discards
  its result. It is cited as the reference for the vendored physics, so it should
  be corrected or retired, not left to mislead.

## 3. Literature gate — OPEN

**Blocking.** Questions:

1. Which closure for which particle class, and on what criterion — Reynolds
   number, cohesiveness, shape? What does a per-element dispatch actually need to
   choose between?
2. Is a *fixed* fractal dimension defensible, or does `Df` have to vary with floc
   size (which changes effective density, and hence settling, non-trivially)?
3. What is the right drag law for genuinely irregular solids (debris, paint
   chips, plastic fragments), beyond the Corey-factor correction?
4. Does particle density evolve over a multi-year run (biofouling), and does that
   matter for a pollutant carrier?
5. What does "offline setup" need to expose so a user can build a defensible
   particle-class spectrum without running the model?

Candidate papers — **all unvetted**:

*Solid grains and shape*
- Stokes; Dietrich (1982); Ferguson & Church (2004)
- Bagheri & Bonadonna (2016) *(implemented)*; Corey (1949)
- Dioguardi et al. (2018), drag for irregular particles over a wide Re range

*Cohesive flocs*
- Maggi (2013) *(implemented; data on disk)*; Kranenburg (1994)
- Khelifa & Hill (2006), variable fractal dimension — **the named gap**
- Winterwerp (1998, 2002); Strom & Keyvani (2011); Manning & Dyer (2007)

*Non-mineral carriers*
- Waldschläger & Schüttrumpf (2019); Khatmullina & Isachenko (2017)
- Kooi et al. (2017), biofouling and vertical transport

## 4. Design

*Stub — pending §3.*

## 5. Implementation steps

*Stub — pending §3.*

## 6. Validation

Existing: Maggi (2013) measurements, the Stokes low-Re limit, BB16 monotonicity in
sphericity — all in `tests/models/test_settling.py`. To add: per-element dispatch
consistency (one run of mixed classes must equal the union of single-class runs),
and, once `settling_dynamic` is on, that a static run with constant state
reproduces the dynamic one exactly.

## 7. Open questions & decisions log

- Does `settling_dynamic` recomputing every step cost enough to matter at 162k
  particles? It was left off partly for that reason and never measured.
- `sed_class` currently conflates cohesiveness with a seeding-class index; M3
  should split it (see ARCHITECTURE.md, element state).

## 8. Status

**Planned. Literature gate open — do not start §4.**
