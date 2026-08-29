# SPPMDRIFT — suspended polluted particulate matter drift

A fully 3D Lagrangian model that **decouples the sediment bed from the water
column** and tracks the **centre of mass of a pollutant** as it moves between
particulate and dissolved form — rather than tracking an inert sediment tracer.

It merges what `SedimentDrift` knows about the bed (wave–current bed stress,
dynamic erosion thresholds, physics-based settling) with what `ChemicalDrift`
knows about chemistry (partitioning, degradation, volatilization). Neither model
can do this today, and each is strong exactly where the other is weak.

> **Scope.** This is a *separate model-development project*. It grew out of the
> DDT dispersal study in `ddt_dump`, but it is not that study and must not be
> driven by it. See [CLAUDE.md](CLAUDE.md).

## Status board

| # | Module | State | Gate |
|---|--------|-------|------|
| 0 | [Foundations](M0_foundations.md) — port, tests, fixtures | **done** 2026-08-29 | — |
| 1 | [Resuspension validation](M1_resuspension_validation.md) | **in progress** | literature gate cleared |
| 2 | [Burial & resuspension thresholds](M2_burial_and_thresholds.md) | planned | **literature gate open** |
| 3 | [Settling](M3_settling.md) | planned | **literature gate open** |
| 4 | [Chemistry & the dissolved fraction](M4_chemistry_dissolved.md) | planned | **literature gate open** |
| 5 | [Dynamic particle size](M5_dynamic_size.md) | proposed | **go/no-go pending** |

Modules are ordered by difficulty, lowest first. They are *not* strictly
sequential — M2 and M3 are largely independent — but M4 depends on M2/M3 for the
bed and carrier physics, and **M5 branch B strictly depends on M4** (see below).

## The working agreement

Each module document follows one skeleton:

```
1. Goal & scope          what "done" means, and what is explicitly out of scope
2. Current state         what already exists, with file:line
3. Literature gate       BLOCKING. Compile and review the papers, jointly.
4. Design                stub until §3 clears
5. Implementation steps  stub until §3 clears
6. Validation            how we know it works
7. Open questions & decisions log
8. Status
```

**No module proceeds past §3 until the paper compilation has been done
together.** The candidate reference lists in §3 and in [REFERENCES.md](REFERENCES.md)
are *proposals*, marked `unvetted`, to be pruned and extended — not a reading
list to be implemented from.

Two further standing rules, from how this project is run:

- **Validate physics against published idealized cases, never against the parent
  study's observational deposition data.** An observed deposition pattern
  confounds transport physics with source history and cannot isolate a physics
  defect. See M1 §3.
- **One module at a time.** The whole plan set is drafted up front so the
  architecture is coherent, but only one module is executed at a time.

## Architecture in one paragraph

`SPPMDrift(OceanDrift)` in `opendrift/models/sppmdrift.py`, with the physics
factored into a dependency-free `opendrift/models/sppm/` package that
`SedimentDrift` also imports. `SedimentDrift` stays behaviour-identical, so the
DDT production runs already in flight remain reproducible. Details, including the
separation-of-concerns contract every module must respect, are in
[ARCHITECTURE.md](ARCHITECTURE.md).

## Where things live

| | |
|---|---|
| Model code | `opendrift/models/` (branch `sppmdrift`) |
| Physics tests | `tests/models/test_{bblm_sg2000,floc_strength,rouse,settling,sedimentdrift_bbl}.py` |
| Test fixtures | `tests/test_data/sppm/` (+ `PROVENANCE.md`) |
| Validation report | `docs/sppmdrift/validation/REPORT.md` |
| Papers | `ddt_dump/resuspension_mechanics/`, `ddt_dump/bblm_sg2000/` (external) |

## Running the physics tests

```
/home/smullersoares/anaconda3/envs/opendrift_dev/bin/python -m pytest \
    tests/models/test_bblm_sg2000.py tests/models/test_floc_strength.py \
    tests/models/test_rouse.py tests/models/test_settling.py \
    tests/models/test_sedimentdrift_bbl.py -o addopts="" -q
```

`opendrift_dev` is the environment with this repository installed editable;
`-o addopts=""` skips the project's `--benchmark-disable --doctest-modules`
defaults, which need plugins that are not installed here.
