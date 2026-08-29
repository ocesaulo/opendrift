# SPPMDRIFT architecture

## Target layout

```
opendrift/models/sppm/            pure-NumPy physics; imports nothing from opendrift
    settling.py                   moved  - terminal-velocity closures
    bblm_sg2000.py                moved  - Styles & Glenn (2000) wave-current BBL
    floc_strength.py              moved  - Shields / fractal-floc critical stress
    rouse.py                      moved  - near-bed Rouse height sampling
    bed.py                        NEW M2 - burial, consolidation, exposure
    partition.py                  NEW M4 - Kd/KOC equilibrium + kinetic sorption
    exchange.py                   NEW M4 - colocation-kernel mass transfer
    aggregation.py                NEW M5 - floc kinetics / population balance

opendrift/models/sedimentdrift.py        behaviour unchanged; imports from sppm/
opendrift/models/sppmdrift.py            NEW - class SPPMDrift(OceanDrift)
```

The four existing modules are pure functions over NumPy arrays already, with no
OpenDrift imports. That is the property worth preserving: it is what makes them
testable without a reader, a run, or a forcing file, and it is why the M1
benchmark suite can be fast enough for CI.

### Two consequences of the move, to handle deliberately

1. **Provenance breakage.** `ddt_dump/codes/ddt_sim/run.py:93` records an md5 of
   `settling.py`, `floc_strength.py` and `rouse.py` **by path**, into every run's
   JSON sidecar. Moving them silently breaks that recording. Either keep
   re-export shims at the old paths for one release, or update `run.py` in the
   same change. Do not do the move without deciding which.
2. **Back-compat imports.** `from opendrift.models import settling` appears in
   the DDT project's code. Shims cost two lines each and remove the whole class of
   breakage; prefer them.

## The separation-of-concerns contract

This is the rule that keeps the model from tangling as chemistry is added. It
generalises the split the SG2000 work already established ("bed stress is a
property of the environment and the bed, never of the particle you happen to be
tracking"), which is why `bed_median_grain_size`, `bed_sediment_density` and
`bed_porosity` are config rather than element variables.

| # | Quantity | Is a function of | Never of |
|---|----------|------------------|----------|
| 1 | Bed shear stress | environment (currents, waves, T/S, depth) + bed properties | the tracked particle |
| 2 | Erosion threshold, settling velocity | the particle: size, density, cohesiveness, shape, floc state | the bed stress |
| 3 | Bed state (burial depth, consolidation) | the particle's own history + the bed environment | other particles (no Eulerian bed) |
| 4 | Partitioning | compound x carrier x ambient chemistry (SPM, DOC, f_OC, pH, T) | position |
| 5 | Transport | the phase mass distribution | — *(new in M4)* |

Rule 5 is the one that changes the model's character. Today an element has one
velocity. Once an element carries pollutant mass in more than one phase, "the
velocity of the element" stops being well defined, because the particulate
fraction settles and the dissolved fraction does not. How to resolve that is
M4's central open question and is deliberately unresolved here.

## Element state, staged

Modules 1–3 should reserve the slots M4/M5 need rather than reshaping the element
class twice.

**Present today** (`SedimentElement`): `settled`, `moving`, `terminal_velocity`,
`grain_diameter`, `d0`, `rho_s`, `corey_shape_factor`, `tau_crit`, `sed_class`,
`fractal_dim`, `phi0`, `time_since_settled`, `C_l`, `use_stokes`, `vel_set`,
`times_resuspended`, `beached`, `latest_resuspension_height`, `viscosity_molecular`.

**Added by module**

| Module | Variable | Purpose |
|--------|----------|---------|
| M1 | `mass` | so a particle cloud can be compared with an Eulerian concentration in kg m^-3 (M4 needs it regardless) |
| M2 | `burial_depth` | depth below the sediment–water interface; gates resuspension |
| M2 | `bed_residence_time` | distinct from `time_since_settled` if consolidation becomes overburden-driven |
| M3 | `cohesive` | explicit flag, replacing the inferred `sed_class` / 63 µm cutoff |
| M3 | `settling_model` | per-element closure code, so one run can host flocs and solid debris |
| M4 | `mass_phase[n_phase]` | the pollutant mass vector — replaces a scalar mass plus an integer phase |
| M4 | *(TBD)* | whatever the dissolved-representation decision requires |
| M5 | — | reuses `grain_diameter`, `d0`, `fractal_dim`; they simply stop being constant |

**Naming note.** `sed_class` currently means two things — a cohesiveness flag and
a seeding-class index. M3 should split them; `cohesive` above is that split.

## Why a mass vector, not a specie index (M4 preview)

`ChemicalDrift` gives each element an integer `specie` and moves it wholesale
between phases by drawing a random number
(`chemicaldrift.py:1415 update_partitioning`), including a Python `for` loop over
every transforming element. The proposed alternative is that each element carries
`mass_phase[n_phase]` and the *same* rate matrix that `init_transfer_rates`
already builds integrates it deterministically.

This is not a stylistic preference. It (a) removes Monte-Carlo partitioning
noise, so particles buy transport statistics instead of partitioning statistics;
(b) removes an O(N) Python loop that will not survive the 162k-particle runs this
lineage of model already does; and (c) is literally "the centre of mass of the
pollutant", which is the project's stated goal. The formulation has a literature
anchor in mass-transfer particle tracking. See M4 §3.

## Relationship to `SedimentDrift` and `ChemicalDrift`

`SedimentDrift` and `ChemicalDrift` are strong in complementary places, and
SPPMDRIFT is the union:

| | `SedimentDrift` | `ChemicalDrift` |
|---|---|---|
| Bed stress | SG2000 wave–current BBL, guard by grid resolution | none |
| Erosion threshold | dynamic Shields / fractal floc + consolidation | a single depth-averaged current speed (`resuspension_critvel`, 0.01 m/s) |
| Resuspension height | Rouse draw, timestep-consistent | fixed depth + Gaussian noise |
| Settling | 5 closures, shape and fractal aware | one hard-coded Sundby Stokes law |
| Burial | none | a fixed matrix rate, built once at seed time |
| Chemistry | none | partitioning, degradation, volatilization, 12 species |

The symmetry points are clean — `bottom_interaction` (`sedimentdrift.py:753` /
`chemicaldrift.py:1581`), `resuspension` (`:838` / `:1607`),
`update_terminal_velocity` (`:685` / `:1114`) — but `ChemicalDrift`'s element
class has only `diameter` and `density`, and its `required_variables` lack the
wave fields SG2000 needs. The port is therefore not a copy; it is a merge.
