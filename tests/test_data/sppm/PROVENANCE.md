# SPPM test fixtures — provenance

Small reference data for the sediment/BBL/settling physics added on the
`sppmdrift` branch (`opendrift/models/{bblm_sg2000,floc_strength,rouse,settling}.py`).

## `test2.mat` (208 B)

Input table for the Styles & Glenn (2000) combined wave–current bottom boundary
layer user test case, distributed with the authors' MATLAB implementation
(`UserTestCase.pdf`, Table 2). Columns used: `Ub`, `Ab`, `Ur`, `zr`, `deg` (cgs).

Expected outputs (cm s^-1 and cm), hard-coded in
`tests/models/test_bblm_sg2000.py`:

| case  | ustarcw | ustarwm | ustarc | znot | z1  |
|-------|---------|---------|--------|------|-----|
| SW-SC | 12.9    | 11.2    | 6.5    | 1.3  | 3.0 |
| WW-SC | 6.7     | 3.0     | 6.0    | 0.1  | 1.4 |
| SW-WC | 17.1    | 17.0    | 1.1    | 0.6  | 3.1 |

Reference: Styles, R. and Glenn, S. M. (2000), Modeling stratified wave and
current bottom boundary layers on the continental shelf, *J. Geophys. Res.*
105(C10), 24119–24139.

## `matlab_sg2000/` (7 files, ~17 kB)

The authors' original MATLAB implementation, from which
`opendrift/models/bblm_sg2000.py` was ported. Kept so the port can be re-derived
and audited:

- `bblm02.m`    — driver
- `bstress2.m`  — bed-stress kernel
- `phi2_1.m`    — Kelvin-function velocity-defect function
- `pwave.m`     — pure-wave limiting solution
- `shldc.m`     — Shields curve
- `neutsed_fun.m`, `CurSedPro.m` — neutral sediment / current profiles
  (**not yet ported**; retained for the deferred profile work)

## `df_particles.csv` (64 kB)

Measured settling velocities of solid mineral grains compiled by Maggi (2013),
used as the observational check on the settling closures
(`tests/models/test_settling.py`). Columns used: `L_microns`, `v_mm_s`, `rho_s`,
`rho_w_fit`, `delta` (only `delta == 3`, i.e. solid non-fractal grains).

Reference: Maggi, F. (2013), The settling velocity of mineral, biomineral, and
biological particles and aggregates in water, *J. Geophys. Res. Oceans* 118,
2118–2132, doi:10.1002/jgrc.20086.

Upstream copies of all of the above live outside this repository, in the
`ddt_dump` project (`bblm_sg2000/`, `data/maggi_test_data/`).
