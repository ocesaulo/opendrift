# M0 — Foundations

**Status: done, 2026-08-29.** Recorded here because it changed what the branch
contains and corrected two live defects.

## 1. Goal & scope

Make the fork self-contained and trustworthy before any new physics is planned:
one source of truth for the code, tests that ship with the branch, and no
docstrings pointing at files that do not exist.

Out of scope: any change to `SedimentDrift`'s numerical behaviour.

## 2. What was found

The premise going in was that functionality lived only in the local executable
install and needed porting. **It did not.** A recursive diff of
`lana_opendrift/.../opendrift/` against the fork returned exactly one differing
file, `models/basemodel/environment.py`, and the whole difference was nine lines
of *commented-out* debug prints. `ddt_dump/bblm_sg2000/sg2000.py` was
byte-identical to `opendrift/models/bblm_sg2000.py`.

The real gap was everything around the code:

- the branch shipped 2164 lines of new physics that `tests/` never touched;
- `test2.mat` and the MATLAB sources that `bblm_sg2000.py`'s docstring said were
  "in this folder" were not in the repository at all;
- `test_sg2000.py` imported a *local copy* of the solver, so the Table 2 result
  did not validate shipped code;
- 12 commits of physics work were unpushed.

## 3. What was done

| | |
|---|---|
| Pushed | 12 commits to `origin/ll_sediment_improvements`; branched `sppmdrift` from it |
| Fixtures | `tests/test_data/sppm/` — `test2.mat`, `matlab_sg2000/`, `df_particles.csv`, `PROVENANCE.md` |
| Tests | 30 tests across `test_{bblm_sg2000,floc_strength,rouse,settling,sedimentdrift_bbl}.py` |
| Fixed | legacy `calc_bottom_stress` ignored `vertical_mixing:bottom_drag_coefficient` |
| Fixed | `calc_upward_resuspension_velocity` read non-existent `elements.E_0` / `elements.porosity` |
| Added | configs `vertical_mixing:erosion_rate_constant`, `vertical_mixing:bed_porosity` |
| Fixed | `examples/example_sediments_resuspension.py` set a config no longer read |
| Corrected | docstrings citing files outside the repository |

The BBL guard test previously needed a real MITgcm run purely to populate wave
fields; it now fabricates the environment and needs no forcing file. That
fabrication helper (`tests/models/sppm_helpers.py`) is the seed of M1's harness.

Not ported: the nine commented-out debug lines in `basemodel/environment.py` —
no functionality.

## 4. Verification

- 30 new tests pass.
- Full `tests/models` suite: 165 passed, 1 failed, 6 errors — **identical to the
  pristine baseline** with the changes stashed, so no regression. (The failure is
  a GDAL geometry API incompatibility in `test_seed_shapefile`; the six errors are
  the ADIOS oil database, which needs network. Both pre-date this work.)
- The repaired example runs and resuspension fires with the seeded `tau_crit`.

## 5. Left undone

- `pytest` had to be installed into `opendrift_dev`; it is still absent from
  `lana_opendrift`, so the suite cannot be run there.
- The `sppm/` package reorganisation in `ARCHITECTURE.md` is **not** done. It
  carries the `ddt_sim/run.py:93` provenance-breakage caveat and should be a
  deliberate, separate change.
- `tests/models/test_{integration,perf}.py` from `ddt_dump/bblm_sg2000/` were not
  ported: they need real MITgcm forcing. M1's harness should replace them.
