# M1 validation notebooks

One notebook per validation case. They are the *illustrated* form of the test
suite: each notebook sets a case up, plots the forcing and the initial condition,
runs it, plots the result against the analytic or published answer, and states a
verdict.

| Notebook | Case | Runtime | Verdict |
|---|---|---|---|
| `00_overview.ipynb` | index + scoreboard | 2 s | — |
| `01_settling_column.ipynb` | settling in isolation | ~45 s | pass |
| `02_well_mixed_condition.ipynb` | the vertical random walk | ~60 s | pass |
| `03_rouse_profile.ipynb` | settling + mixing + lift-off height | ~3 min | pass, documented bias |
| `04_S3_mixed_bed_closure.ipynb` | Sherwood Fig. 2 mixed-bed threshold | ~20 s | **fail (F1)** |
| `05_S1_double_resuspension.ipynb` | Sherwood §3.2.1 double resuspension | ~5 min | 4 pass, **2 fail (F2, F3)** |

Start with `00_overview.ipynb`.

## Regenerating

The notebooks are **generated**, not hand-edited — `build_notebooks.py` is the
source of truth, so a notebook and its test cannot drift apart: both import the
same harness, `tests/models/sppm_harness.py`.

```
cd docs/sppmdrift/validation/notebooks
$ENV/bin/python build_notebooks.py                 # write the .ipynb files
$ENV/bin/python build_notebooks.py --execute       # ... and run them
$ENV/bin/python build_notebooks.py --only 05_ --execute   # just one
```

where `$ENV` is `~/anaconda3/envs/opendrift_dev`, the environment with this
repository installed editable. Edits made in Jupyter will be **overwritten** by
the next build; change `build_notebooks.py` instead.

Execution goes through `nbclient` rather than `jupyter nbconvert --execute`,
because this machine's Jupyter configuration registers a preprocessor from
`jupyter_contrib_nbextensions` that is not installed, so nbconvert fails at
start-up. This mirrors `ddt_dump/codes/run_notebook.py`.

## Relationship to the tests

Each notebook corresponds to tests in
`tests/models/test_sppm_analytic.py` and `tests/models/test_sppm_roms_cases.py`.
The tests are the gate — they run in CI in about 12 seconds with the whole-run
cases skipped, or 8 minutes with `--run-slow`. The notebooks are the explanation.

Narrative write-up: [`../REPORT.md`](../REPORT.md).
