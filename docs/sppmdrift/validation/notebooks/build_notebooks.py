#!/usr/bin/env python
"""Build (and optionally execute) the SPPMDRIFT validation notebooks.

    python build_notebooks.py            # write the .ipynb files
    python build_notebooks.py --execute  # ... and run them, storing outputs

One notebook per validation case. They import the same harness the test suite
uses (`tests/models/sppm_harness.py`), so a notebook and its test cannot drift
apart: the notebook shows what the assertion is asserting.

Execution goes through nbclient rather than `jupyter nbconvert --execute`,
because this machine's Jupyter configuration registers a preprocessor from
`jupyter_contrib_nbextensions` which is not installed here, so nbconvert fails at
start-up. This mirrors ddt_dump/codes/run_notebook.py.
"""
import argparse
import pathlib
import sys
import time

import nbformat as nbf

HERE = pathlib.Path(__file__).resolve().parent
KERNEL = "python3"          # resolves to the env running this script

HEADER = '''\
# {TITLE}
import sys, pathlib, warnings
warnings.filterwarnings("ignore")

# locate the repository root (the directory containing `opendrift/`)
_p = pathlib.Path.cwd().resolve()
while not (_p / "opendrift" / "__init__.py").exists() and _p != _p.parent:
    _p = _p.parent
REPO = _p
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline

plt.rcParams.update({
    "figure.figsize": (9, 4.2), "figure.dpi": 110,
    "axes.grid": True, "grid.alpha": 0.3,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
})
CLR = ["#0b6fa4", "#e07b39", "#3f9b52", "#b4434e", "#7c5aa6"]
print("repo:", REPO)
'''


def nb(cells):
    n = nbf.v4.new_notebook(cells=cells)
    n.metadata.update({
        "kernelspec": {"display_name": "Python 3", "language": "python",
                       "name": "python3"},
        "language_info": {"name": "python"},
    })
    return n


def md(text):
    return nbf.v4.new_markdown_cell(text.strip("\n"))


def code(text):
    return nbf.v4.new_code_cell(text.strip("\n"))


def header(title):
    return code(HEADER.replace("{TITLE}", title))


BUILDERS = {}


def builder(filename):
    def deco(fn):
        BUILDERS[filename] = fn
        return fn
    return deco


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--execute", action="store_true")
    ap.add_argument("--only", default=None, help="substring of the filename")
    ap.add_argument("--timeout", type=int, default=3600)
    args = ap.parse_args()

    targets = {k: v for k, v in BUILDERS.items()
               if args.only is None or args.only in k}
    for name, fn in sorted(targets.items()):
        path = HERE / name
        nbf.write(fn(), path)
        print("wrote", path.name)

    if not args.execute:
        return
    from nbclient import NotebookClient
    for name in sorted(targets):
        path = HERE / name
        print("executing", name, "...", end=" ", flush=True)
        book = nbf.read(path, as_version=4)
        client = NotebookClient(book, timeout=args.timeout, kernel_name=KERNEL,
                                allow_errors=False,
                                resources={"metadata": {"path": str(HERE)}})
        t0 = time.time()
        try:
            client.execute()
        finally:
            nbf.write(book, path)
        print("ok (%.0f s)" % (time.time() - t0), flush=True)


# =========================================================================
@builder("01_settling_column.ipynb")
def _settling_column():
    return nb([
md("""
# Tier A — the settling column

**What this tests.** Settling, in isolation. Turbulent mixing is switched off and
the critical shear stress is set impossibly high, so nothing is resuspended: the
only thing happening is that particles fall at a constant terminal velocity and
stop when they reach the bed.

**The analytic answer.** A cloud spread uniformly through a column of depth $h$,
all settling at the same $w_s$, deposits *linearly* in time. After a time $t$,
everything that started within $w_s t$ of the bed has landed:

$$ f_{\\rm settled}(t) = \\min\\!\\left(1,\\ \\frac{w_s t}{h}\\right) $$

**Why it matters.** This is the cheapest possible check that terminal velocity is
applied with the right sign and magnitude, and that `bottom_interaction` retires
elements exactly when they reach the seafloor — no earlier, no later.

Test: `tests/models/test_sppm_analytic.py::test_settling_column_deposition_rate`
"""),
header("Tier A - settling column"),
md("""
## Setup

A 20 m column, no current, **no turbulent mixing**. 4000 elements are spread
uniformly through the water column and settle at 2 mm/s, so the column is swept
clean in $h/w_s = 10^4\\,$s $\\approx 2.8\\,$h. We run for 2 h, which should leave
about 72 % deposited.
"""),
code("""
import datetime, tempfile, os
import xarray as xr
from tests.models.sppm_harness import column_dataset, column_model, seed_column

H, W_S, HOURS, N = 20.0, 2.0e-3, 2.0, 4000

ds = column_dataset(h=H, u=0.0, u_star=0.02, hours=HOURS + 1)
o, reader = column_model(
    ds,
    drift__vertical_mixing=False,                 # settling in isolation
    vertical_mixing__settling_model='prescribed',
    vertical_mixing__tau_crit_mode='constant')
z0 = seed_column(o, reader, N, H, W_S, tau_crit=1e6)  # tau_crit huge -> no resuspension

print(f"h = {H} m,  w_s = {W_S*1e3:.1f} mm/s,  sweep time h/w_s = {H/W_S/3600:.2f} h")
print(f"{N} elements seeded between z = {z0.min():.2f} and {z0.max():.2f} m")
"""),
md("""
### Initial condition

Uniform in depth — this is what makes the deposition linear in time. Any
clustering here would bend the expected curve.
"""),
code("""
fig, ax = plt.subplots(figsize=(5.2, 4.2))
ax.hist(z0, bins=40, orientation='horizontal', color=CLR[0], alpha=.85)
ax.axhline(-H, color='k', lw=2)
ax.text(0.98, -H + 0.4, 'seafloor', transform=ax.get_yaxis_transform(),
        ha='right', va='bottom', fontsize=9)
ax.set_xlabel('elements per bin'); ax.set_ylabel('z  [m]')
ax.set_title('Initial distribution: uniform through the column')
plt.show()
"""),
md("""
## Run

Output every 6 minutes so the deposition curve is well resolved.
"""),
code("""
out = os.path.join(tempfile.mkdtemp(), 'settling_column.nc')
o.run(time_step=60, time_step_output=360,
      duration=datetime.timedelta(hours=HOURS), outfile=out,
      export_variables=['z', 'settled'], stop_on_error=True)

with xr.open_dataset(out) as d:
    settled = d['settled'].values          # (trajectory, time)
    z_hist = d['z'].values
hours = np.arange(settled.shape[1]) * 360 / 3600
frac = np.nanmean(settled == 1, axis=0)
print(f"final settled fraction: {frac[-1]:.3f}   expected {min(1, W_S*HOURS*3600/H):.3f}")
"""),
md("""
## Result

The modelled curve should lie on the analytic line until the column empties.
"""),
code("""
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))

# left: a sample of trajectories, all with the same slope
for i in np.linspace(0, len(z_hist) - 1, 25).astype(int):
    a1.plot(hours, z_hist[i], color=CLR[0], lw=.7, alpha=.5)
a1.axhline(-H, color='k', lw=2)
a1.plot(hours, -H + W_S * hours * 3600, color=CLR[3], lw=2, ls='--',
        label=f'$w_s$ = {W_S*1e3:.0f} mm/s')
a1.set_xlabel('time [h]'); a1.set_ylabel('z [m]'); a1.set_ylim(-H - 1, 0.5)
a1.set_title('Trajectories: constant fall speed'); a1.legend(loc='lower left')

# right: deposition curve against the analytic result
analytic = np.minimum(1.0, W_S * hours * 3600 / H)
a2.plot(hours, analytic, color='k', lw=2.5, ls='--', label=r'analytic  $w_s t / h$')
a2.plot(hours, frac, color=CLR[0], lw=2, label='model')
a2.fill_between(hours, analytic - .05, analytic + .05, color='k', alpha=.10,
                label='tolerance $\\pm$0.05')
a2.set_xlabel('time [h]'); a2.set_ylabel('fraction settled'); a2.set_ylim(0, 1)
a2.set_title('Deposition is linear in time'); a2.legend(loc='upper left')
a2.annotate(f'final: {frac[-1]:.3f}\\nexpected: {analytic[-1]:.3f}',
            xy=(hours[-1], frac[-1]), xytext=(-110, -40),
            textcoords='offset points', fontsize=9,
            arrowprops=dict(arrowstyle='->', color='0.3'),
            bbox=dict(boxstyle='round,pad=0.4', fc='white', ec='0.7'))
plt.tight_layout(); plt.show()
"""),
md("""
## Verdict

**PASS.** The deposition curve follows $w_s t/h$ within the 0.05 tolerance
throughout. Terminal velocity is applied correctly and elements are retired at
the seafloor exactly when they arrive.
"""),
    ])


# =========================================================================
@builder("02_well_mixed_condition.ipynb")
def _well_mixed():
    return nb([
md("""
# Tier A — the well-mixed condition

**What this tests.** The vertical random walk, in isolation. Settling is switched
off ($w_s = 0$), so elements are neutrally buoyant tracers.

**The requirement.** This is the standard consistency condition on a vertical
random-walk scheme (Thomson 1987; Visser 1997; Ross & Sharples 2004):

> An initially uniform distribution of neutrally buoyant particles must **stay**
> uniform, even when the diffusivity varies by orders of magnitude over the water
> column.

**The classic failure.** A naive random walk that uses the local $K$ without the
$\\partial K/\\partial z$ correction drives particles *into* regions of low
diffusivity, where they then take smaller steps and become stuck. With a
parabolic $K(z)$ — small at both the bed and the surface — that failure shows up
as spurious accumulation at **both** boundaries.

**Which part of the column we can test.** The bed is *not* a reflecting boundary
here: this model settles elements onto it and re-releases them at a fixed
reference height, which is a real boundary treatment but not the one this
condition is posed for. So the bottom tenth is excluded.

The **surface is** a genuine reflecting boundary, so the upper column is fair
game — and that is where the test has teeth, because $K$ collapses by four orders
of magnitude as it approaches the surface. Restricting to the middle of the
column, as a first attempt did, quietly throws away almost all the variation in
$K$ and makes the test far weaker than it looks.

Test: `tests/models/test_sppm_analytic.py::test_well_mixed_condition`
"""),
header("Tier A - well-mixed condition"),
md("""
## Setup

The diffusivity profile is the log-layer parabola

$$ K(z_b) = \\kappa\\, u_* \\, z_b \\left(1 - \\frac{z_b}{h}\\right) $$

with $z_b$ the height above the bed. It vanishes at both boundaries and peaks at
mid-depth — deliberately hostile to a naive scheme.
"""),
code("""
import datetime
from tests.models.sppm_harness import (column_dataset, column_model, seed_column,
                                       height_above_bed, parabolic_diffusivity, KAPPA)

H, U_STAR, HOURS, N = 20.0, 0.02, 12, 6000

zb = np.linspace(0, H, 400)
K = parabolic_diffusivity(H - zb, H, U_STAR)

fig, ax = plt.subplots(figsize=(5.2, 4.2))
ax.plot(K * 1e4, zb, color=CLR[0], lw=2)
ax.axhspan(0, 0.1 * H, color='0.85', zorder=0)
ax.axhspan(0.9 * H, H, color='0.85', zorder=0)
ax.text(K.max() * 1e4 * .55, 0.05 * H, 'excluded', fontsize=8, va='center')
ax.text(K.max() * 1e4 * .55, 0.95 * H, 'excluded', fontsize=8, va='center')
ax.set_xlabel(r'$K$  [$10^{-4}$ m$^2$ s$^{-1}$]'); ax.set_ylabel(r'height above bed $z_b$ [m]')
ax.set_title(r'Parabolic diffusivity: $K\\to0$ at both boundaries')
ax.annotate(r'$K_{max}=\\kappa u_* h/4$', xy=(K.max()*1e4, H/2), xytext=(-95, 22),
            textcoords='offset points', fontsize=9,
            arrowprops=dict(arrowstyle='->', color='0.3'))
plt.tight_layout(); plt.show()
interior = K[(zb > 0.1 * H) & (zb < 0.9 * H)]
print(f"K over the whole column : {K.min():.2e} to {K.max():.2e} m^2/s")
print(f"K over the tested interior: {interior.min():.2e} to {interior.max():.2e} m^2/s"
      f"  ({interior.max()/interior.min():.0f}x range)")
"""),
md("""
## Run

6000 neutrally buoyant elements, initially uniform, for 12 hours.
"""),
code("""
ds = column_dataset(h=H, u=0.0, u_star=U_STAR, hours=HOURS + 1)
o, reader = column_model(
    ds,
    vertical_mixing__settling_model='prescribed',
    vertical_mixing__tau_crit_mode='constant',
    vertical_mixing__resuspension_height_mode='reference',
    vertical_mixing__resuspension_reference_height=0.01)
z_seed = seed_column(o, reader, N, H, w_s=0.0, tau_crit=1e-9)
zb0 = H + z_seed

o.run(time_step=120, time_step_output=3600,
      duration=datetime.timedelta(hours=HOURS), stop_on_error=True)
zb1 = height_above_bed(o, H)
print(f"{len(zb1)} elements still in the domain after {HOURS} h")
"""),
md("""
## Result

If the scheme were biased, the final histogram would sag in the middle and pile
up near the boundaries, tracking the shape of $K$. It should instead stay flat.
"""),
code("""
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2),
                             gridspec_kw={'width_ratios': [1, 1.25]})

bins = np.linspace(0, H, 21)
a1.hist(zb0, bins=bins, orientation='horizontal', histtype='step',
        color='0.4', lw=2, label='initial')
a1.hist(zb1, bins=bins, orientation='horizontal', histtype='stepfilled',
        color=CLR[0], alpha=.55, label=f'after {HOURS} h')
a1.set_xlabel('elements per bin'); a1.set_ylabel(r'$z_b$ [m]')
a1.set_title('Distribution stays uniform'); a1.legend(loc='lower right', fontsize=9)

a2.axhspan(-.20, .20, color=CLR[2], alpha=.13, label='tolerance $\\pm$20 %')
a2.axhline(0, color='k', lw=1)
for (lo_f, hi_f), c, nm in [((0.10, 0.90), CLR[0], 'interior 10-90 %'),
                            ((0.70, 1.00), CLR[1], 'top 30 % ($K\\to0$)')]:
    lo, hi = lo_f * H, hi_f * H
    counts, edges = np.histogram(zb1[(zb1 > lo) & (zb1 < hi)], bins=8, range=(lo, hi))
    centres = 0.5 * (edges[1:] + edges[:-1])
    dev = (counts - counts.mean()) / counts.mean()
    Kr = parabolic_diffusivity(H - np.linspace(lo, hi, 60), H, U_STAR)
    a2.plot(centres, dev, 'o-', color=c, lw=1.8, ms=5,
            label=f'{nm}:  K varies {Kr.max()/max(Kr.min(),1e-12):.0f}x,'
                  f'  max dev {100*np.abs(dev).max():.0f} %')
a2.set_xlabel(r'$z_b$ [m]'); a2.set_ylabel('deviation from uniform')
a2.set_title('No structure in either region')
a2.legend(fontsize=8, loc='lower left'); a2.set_ylim(-.32, .32)
plt.tight_layout(); plt.show()
"""),
md("""
## Verdict

**PASS**, and on the demanding version of the test.

The interior check is the weak one: excluding the outer tenths leaves only about
a **3×** spread in $K$, so passing it says little. The upper-column check is the
real one — across the top 30 % the diffusivity falls by more than **four orders of
magnitude**, and that is exactly where a scheme missing the $\partial K/\partial z$
correction would pile elements up against the surface. The distribution stays
uniform to within ~10 %, comfortably inside tolerance.

The random walk is therefore consistent with a strongly depth-varying
diffusivity, and the mixing underneath the Rouse case in notebook 03 can be
trusted.
"""),
    ])


# =========================================================================
@builder("03_rouse_profile.ipynb")
def _rouse():
    return nb([
md("""
# Tier A — the Rouse equilibrium profile

**What this tests.** Settling, turbulent mixing and the resuspension release
height, *jointly*. This is the most diagnostic case in Tier A: it is the only one
where all three have to be right at once for the answer to come out.

**The analytic answer.** With a parabolic diffusivity and a constant settling
velocity, the steady state of the vertical advection–diffusion balance is the
Rouse (1937) distribution,

$$ \\frac{C(z_b)}{C(a)} = \\left[\\frac{h-z_b}{z_b}\\cdot\\frac{a}{h-a}\\right]^{P},
\\qquad P = \\frac{w_s}{\\kappa u_*} $$

$P$ is the Rouse number: small $P$ means a well-mixed suspension, large $P$ means
material hugging the bed. Recovering $P$ from the modelled profile is the test.

**Why `reference` release, not the production `turbulent` mode.** The `turbulent`
mode draws lift-off heights *from the Rouse distribution itself*, so testing a
Rouse profile against it would be partly circular. Using a fixed release height
removes that: the profile then has to be built by settling against mixing.

Test: `tests/models/test_sppm_analytic.py::test_rouse_equilibrium_profile`
"""),
header("Tier A - Rouse equilibrium profile"),
md("""
## Setup

A 20 m column with $u_* = 0.02$ m/s. We sweep two Rouse numbers by choosing
$w_s = P \\kappa u_*$; each is a separate 24 h run.
"""),
code("""
import datetime
from tests.models.sppm_harness import (column_dataset, column_model, seed_column,
                                       height_above_bed, rouse_profile, KAPPA)

H, U_STAR, ZA, N, HOURS = 20.0, 0.02, 0.01, 4000, 24
PS = [0.25, 0.50]

for P in PS:
    print(f"P = {P:4.2f}  ->  w_s = {P*KAPPA*U_STAR*1e3:5.2f} mm/s")

zb = np.logspace(np.log10(ZA), np.log10(H), 300)
fig, ax = plt.subplots(figsize=(5.4, 4.2))
for P, c in zip([0.1, 0.25, 0.5, 1.0, 2.0], CLR):
    ax.plot(rouse_profile(zb, P, H, ZA), zb, color=c, lw=2, label=f'P = {P}')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel(r'$C/C(a)$'); ax.set_ylabel(r'height above bed $z_b$ [m]')
ax.set_title('The Rouse family: larger $P$ hugs the bed')
ax.legend(fontsize=8.5)
plt.tight_layout(); plt.show()
"""),
md("""
## Run

Each run seeds a uniform cloud, sets `tau_crit` effectively to zero so material
recirculates continuously, and lets the profile find its equilibrium over 24 h.
"""),
code("""
results = {}
for P in PS:
    w_s = P * KAPPA * U_STAR
    ds = column_dataset(h=H, u=0.3, u_star=U_STAR, hours=HOURS + 2)
    o, reader = column_model(
        ds,
        vertical_mixing__settling_model='prescribed',
        vertical_mixing__tau_crit_mode='constant',
        vertical_mixing__resuspension_height_mode='reference',
        vertical_mixing__resuspension_reference_height=ZA,
        vertical_mixing__bbl_scheme='legacy')
    seed_column(o, reader, N, H, w_s, tau_crit=1e-9)
    o.run(time_step=300, time_step_output=3600,
          duration=datetime.timedelta(hours=HOURS), stop_on_error=True)
    results[P] = dict(zb=height_above_bed(o, H),
                      n_bed=int(np.sum(o.elements.settled == 1)))
    print(f"P={P}: {results[P]['n_bed']} of {N} elements resting on the bed at the end")
"""),
md("""
## Result

Concentration is a count per unit height in log-spaced bins. Elements *resting on
the bed* awaiting their next lift-off are excluded — see the caveat below.
"""),
code("""
fig, axes = plt.subplots(1, len(PS) + 1, figsize=(13, 4.3),
                         gridspec_kw={'width_ratios': [1] * len(PS) + [0.9]})
fits = {}
for ax, P in zip(axes, PS):
    zbv = results[P]['zb']
    edges = np.logspace(np.log10(ZA), np.log10(H), 25)
    counts, _ = np.histogram(zbv, bins=edges)
    centres = np.sqrt(edges[1:] * edges[:-1])
    conc = counts / np.diff(edges)
    ok = counts > 20
    Pfit = np.polyfit(np.log((H - centres[ok]) / centres[ok]), np.log(conc[ok]), 1)[0]
    fits[P] = Pfit

    ref = rouse_profile(centres, P, H, ZA)
    scale = conc[ok][0] / ref[ok][0]
    ax.plot(ref * scale, centres, color='k', lw=2.5, ls='--', label=f'Rouse, P = {P}')
    ax.plot(conc[ok], centres[ok], 'o-', color=CLR[0], ms=4, lw=1.6, label='model')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('concentration [elements / m]'); ax.set_ylabel(r'$z_b$ [m]')
    ax.set_title(f'P = {P}:  fitted {Pfit:.3f}')
    ax.legend(fontsize=8.5, loc='lower left')

ax = axes[-1]
x = np.arange(len(PS))
ax.bar(x - .18, PS, .36, color='0.35', label='true P')
ax.bar(x + .18, [fits[P] for P in PS], .36, color=CLR[0], label='fitted P')
for i, P in enumerate(PS):
    ax.text(i + .18, fits[P] + .01, f'{fits[P]:.3f}', ha='center', fontsize=8.5)
ax.set_xticks(x); ax.set_xticklabels([f'P = {P}' for P in PS])
ax.set_ylabel('Rouse exponent'); ax.set_title('Recovered exponent (biased low)')
ax.legend(fontsize=8.5)
plt.tight_layout(); plt.show()
for P in PS:
    print(f"P = {P}: fitted {fits[P]:.3f}  ({100*(fits[P]-P)/P:+.0f} %)")
"""),
md("""
## Verdict, and a documented bias

**PASS** — the fitted exponent lands within the 0.55–1.15 $\\times P$ tolerance,
and the modelled profile tracks the Rouse shape over more than two decades of
concentration. The settling-versus-mixing balance is right.

**But the estimate is systematically low**, and the reason is worth recording.
Elements resting on the bed between lift-offs are not part of the suspended
profile, so near-bed concentration is under-counted and the profile flattens. A
sweep over how that population is treated:

| $P$ | outer $\\Delta t$ | duration | bed-resting | fitted $P$ |
|---|---|---|---|---|
| 0.25 | 300 s | 24 h | 267 | **0.223** (excluded — what we do) |
| 0.25 | 300 s | 24 h | 267 | 0.582 (counted at $z_a$) |
| 0.25 | 60 s | 24 h | 54 | 0.407 (counted at $z_a$) |
| 0.25 | 60 s | 48 h | 47 | 0.443 (counted at $z_a$) |

Counting them at the release height over-corrects badly; excluding them is much
the better estimator. Shrinking the time step shrinks the bed-resting population
as expected (267 → 54 for a 5× smaller step), which confirms this is a
**discretisation artifact of the lift-off cycle**, not a physical result.

At $P \\gtrsim 1$ so much material sits on the bed that the suspended population
becomes too small to fit, so this benchmark is scoped to $P \\lesssim 0.5$.
"""),
    ])


# =========================================================================
@builder("04_S3_mixed_bed_closure.ipynb")
def _s3():
    return nb([
md("""
# Tier B / S3 — the mixed-bed erosion threshold

**Reference.** Sherwood, C. R. *et al.* (2018), *Cohesive and mixed sediment in
the Regional Ocean Modeling System (ROMS v3.6) implemented in COAWST*,
Geosci. Model Dev. **11**, 1849–1871 — Section 2.6 and Figure 2.

**What this tests.** How the erosion threshold of a *mixed* bed — sand and mud
together — depends on how much mud is present. This is a pure closure comparison:
no simulation is run, so it is instant, and it isolates one equation.

**The physical picture.** A bed with only a few per cent mud behaves like sand:
each grain is held by its own weight, and fines are winnowed out easily. Above
roughly 3–30 % mud the bed starts behaving as a *bulk* cohesive material: the mud
matrix binds everything, so even coarse grains are held by the mud's strength
rather than their own.

**Sherwood's Eq. 6:**

$$ \\tau_{ce} = \\max\\!\\left[\\,P_c\\,\\tau_{cb} + (1-P_c)\\,\\tau_c\\ ,\\ \\tau_c\\,\\right] $$

with $P_c$ the cohesive-behaviour parameter (0 = non-cohesive, 1 = cohesive),
$\\tau_{cb}$ the **bulk** critical stress of the bed — *one number for every size
class* — and $\\tau_c$ the particle's own Shields stress. Note the `max`: making a
bed muddier can never make a grain *easier* to erode.

**What we do instead** (`sedimentdrift.py`, mode `mixed`):

$$ \\tau_{crit} = (1-P_c)\\,\\tau_{\\rm Shields}(d) + P_c\\,\\tau_{\\rm floc}(d, D_f, \\phi) $$

Both terms are per-particle, and $\\tau_{\\rm floc} \\propto d$. This notebook shows
what that difference costs.

Test: `tests/models/test_sppm_roms_cases.py` (the `test_s3_*` family)
"""),
header("Tier B / S3 - mixed-bed erosion threshold"),
md("""
## Setup

The four sediment classes of Sherwood's experiments, and the Figure 2
configuration: bulk bed stress $\\tau_{cb} = 0.1$ Pa, applied bottom stress
$\\tau_b = 0.12$ Pa.
"""),
code("""
from tests.models.test_sppm_roms_cases import (
    CLASS_D, CLASS_TAU_CE, CLASS_WS, TAU_CB, TAU_B,
    sherwood_tau_ce, _shields, _model_mixed_tau)

NAMES = ['clay 4 µm', 'silt 30 µm', 'sand 62.5 µm', 'sand 140 µm']
tau_c = _shields(CLASS_D)

print(f"{'class':>14s} {'d [µm]':>8s} {'w_s [mm/s]':>11s} {'Shields τ_c [Pa]':>17s}")
for n, d, w, t in zip(NAMES, CLASS_D, CLASS_WS, tau_c):
    print(f"{n:>14s} {d*1e6:8.1f} {w*1e3:11.2f} {t:17.4f}")
print(f"\\nbulk bed stress τ_cb = {TAU_CB} Pa;  applied stress τ_b = {TAU_B} Pa")
"""),
md("""
### The cohesive-behaviour ramp

$P_c$ rises from 0 to 1 as the bed's mud fraction increases. Our ramp runs over
5 %→35 % mud (from Jacobs *et al.* 2011 and Yao *et al.* 2022); Sherwood's Figure 2
example completes at 20 %. That difference is minor next to what follows.
"""),
code("""
fmud = np.linspace(0, 0.6, 200)
Pc_ours = np.clip((fmud - 0.05) / (0.35 - 0.05), 0, 1)
Pc_sher = np.clip(fmud / 0.20, 0, 1)

fig, ax = plt.subplots(figsize=(5.6, 3.6))
ax.plot(fmud, Pc_ours, color=CLR[0], lw=2.2, label='ours (5 % → 35 %)')
ax.plot(fmud, Pc_sher, color='k', lw=2, ls='--', label='Sherwood Fig. 2 (→ 20 %)')
ax.set_xlabel('mud fraction  $f_c$'); ax.set_ylabel('$P_c$')
ax.set_title('Cohesive-behaviour parameter')
ax.legend(fontsize=9); plt.tight_layout(); plt.show()
"""),
md("""
## Result — the two formulations side by side

Left: Sherwood Eq. 6. Every class converges toward the *same* bulk bed stress,
floored at its own Shields value, so the curves fan into a narrow band.

Right: ours. Because the cohesive term scales with each particle's **own**
diameter, the classes fan *apart* — a 140 µm grain in a muddy bed is assigned the
fictitious strength of a 140 µm floc.
"""),
code("""
fm = np.linspace(0, 0.6, 25)
ours = np.array([_model_mixed_tau(f, CLASS_D) for f in fm])
sher = np.array([sherwood_tau_ce(np.clip(f / 0.20, 0, 1), TAU_CB, tau_c) for f in fm])

fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.4), sharey=True)
for k, (n, c) in enumerate(zip(NAMES, CLR)):
    a1.plot(fm, sher[:, k], color=c, lw=2, label=n)
    a2.plot(fm, ours[:, k], color=c, lw=2, label=n)
for a, t in ((a1, 'Sherwood Eq. 6 (reference)'), (a2, 'ours: mode "mixed"')):
    a.axhline(TAU_CB, color='0.4', lw=1.2, ls=':')
    a.text(0.60, TAU_CB * 1.06, r'$\\tau_{cb}$', ha='right', fontsize=9, color='0.35')
    a.axhline(TAU_B, color=CLR[3], lw=1.2, ls='-.')
    a.text(0.60, TAU_B * 1.06, r'$\\tau_b$ applied', ha='right', fontsize=9, color=CLR[3])
    a.set_yscale('log'); a.set_xlabel('mud fraction  $f_c$'); a.set_title(t)
a1.set_ylabel(r'effective critical stress  $\\tau_{ce}$  [Pa]')
a1.legend(fontsize=8.5, loc='upper left')
a2.annotate('classes fan apart:\\nheld by their own size,\\nnot by the mud',
            xy=(0.40, ours[-6, 3]), xytext=(-30, -75), textcoords='offset points',
            fontsize=9, ha='center',
            arrowprops=dict(arrowstyle='->', color=CLR[3]),
            bbox=dict(boxstyle='round,pad=0.4', fc='#fff3f0', ec=CLR[3]))
plt.tight_layout(); plt.show()
"""),
md("""
### How large is the error?

In the fully cohesive limit ($P_c = 1$), against Eq. 6:
"""),
code("""
o_lim = _model_mixed_tau(0.5, CLASS_D)
s_lim = sherwood_tau_ce(1.0, TAU_CB, tau_c)
ratio = o_lim / s_lim

fig, ax = plt.subplots(figsize=(6.6, 3.9))
cols = [CLR[2] if 0.8 < r < 1.25 else CLR[3] for r in ratio]
ax.bar(NAMES, ratio, color=cols, alpha=.9)
ax.axhline(1, color='k', lw=1.5)
ax.axhspan(0.8, 1.25, color=CLR[2], alpha=.15, label='within 25 %')
ax.set_yscale('log'); ax.set_ylabel(r'ours / Sherwood Eq. 6')
ax.set_title(r'Error in the fully cohesive limit ($P_c=1$)')
for i, r in enumerate(ratio):
    ax.text(i, r * 1.12, f'{r:.1f}×', ha='center', fontsize=9.5, fontweight='bold')
ax.legend(fontsize=9); plt.tight_layout(); plt.show()

print(f"{'class':>14s} {'Eq.6 [Pa]':>10s} {'ours [Pa]':>10s} {'ratio':>7s}")
for n, s_, o_, r in zip(NAMES, s_lim, o_lim, ratio):
    print(f"{n:>14s} {s_:10.3f} {o_:10.3f} {r:6.1f}×")
"""),
md("""
## Verdict — **FAIL (finding F1)**

Two distinct defects, both visible above:

1. **The cohesive term is the wrong quantity.** It should be the bulk stress of
   the bed, identical for every class. Ours is a floc strength computed from the
   *particle's own* diameter, so coarse grains come out up to **14× too hard to
   erode**.
2. **The `max(·, τ_c)` floor is missing.** Without it, making a bed muddier can
   make the finest class *easier* to erode — the clay bar sits at **0.6×**, below
   one.

**What still works:** at low mud fraction the two agree exactly (both reduce to
the particle's Shields stress), and our thresholds do increase monotonically with
mud content. The structure is right; the cohesive limit is not.

**Scope.** `mixed` is not the production default — `auto` is — so the DDT
production runs are not affected *through this path*. But `auto` uses the same
`tau_crit_floc(d)` for everything below its 63 µm cutoff, which is how a 31 µm
quartz grain in those runs acquires a 0.48 Pa threshold against a legacy value of
0.09 Pa. That is consistent with the near-zero resuspension already seen there.

Fix belongs to **M2/B3**. This is pinned in the test suite as a strict `xfail`,
so a fix cannot land without updating the documentation.
"""),
    ])


# =========================================================================
@builder("05_S1_double_resuspension.ipynb")
def _s1():
    return nb([
md("""
# Tier B / S1 — the double resuspension experiment

**Reference.** Sherwood, C. R. *et al.* (2018), Geosci. Model Dev. **11**,
1849–1871 — Section 3.2.1 and Figure 5.

**The case.** A 20 m column, one-dimensional in the vertical, flat bottom, no floc
dynamics. Two bottom-stress events about 1.5 days apart, lasting 1.5 and 1 days,
the first peaking near 1 Pa. Four sediment classes, all initially in the bed.
Five days total.

**What Sherwood reports, and what we should reproduce:**

1. The finer fractions dominate the suspended load; the water column contains
   *"only a small fraction of the coarsest sand"*.
2. When the stress subsides, *"coarser sediment deposited first, while finer
   material remained suspended"* — building a fining-upward storm layer over a
   coarse lag.
3. *"The finest material (4 µm) remained mostly in suspension after 5 days."*
4. The second, weaker pulse *"only resuspended minimal amounts of the 140 µm
   sand"*.

**The adaptation.** ROMS has a layered bed and reports stratigraphy; we have no
bed layers. So the stratigraphic result — graded bedding, a lag layer — is
compared through **the order in which classes leave suspension**, which is what
builds it. Where a quantity is out of reach we say so rather than approximate it.

Tests: `tests/models/test_sppm_roms_cases.py` (the `test_s1_*` family)
"""),
header("Tier B / S1 - double resuspension experiment"),
md("""
## Setup

### The four sediment classes

Diameters, critical stresses and settling velocities exactly as prescribed in
Section 3.2.1. Note the classes span an **80-fold** range of settling velocity —
that spread is what the deposition-order test exploits.
"""),
code("""
import datetime, tempfile, os
import xarray as xr
from tests.models.sppm_harness import column_dataset, column_model, stress_to_speed
from tests.models.test_sppm_roms_cases import CLASS_D, CLASS_TAU_CE, CLASS_WS, _s1_stress

NAMES = ['clay 4 µm', 'silt 30 µm', 'sand 62.5 µm', 'sand 140 µm']
H, HOURS, N_PER_CLASS = 20.0, 120.0, 500

print(f"{'class':>14s} {'d [µm]':>8s} {'τ_ce [Pa]':>10s} {'w_s [mm/s]':>11s} {'h/w_s [h]':>10s}")
for n, d, t, w in zip(NAMES, CLASS_D, CLASS_TAU_CE, CLASS_WS):
    print(f"{n:>14s} {d*1e6:8.1f} {t:10.2f} {w*1e3:11.1f} {H/w/3600:10.1f}")
"""),
md("""
### The bottom-stress history

Two half-sine events. The forcing is imposed through the current: with the
`legacy` bed-stress scheme and a finely resolved column, the stress reduces to
$\\tau_b = \\rho\\,c_d\\,|u|^2$, so a prescribed stress history is exactly
invertible into a speed.

The dashed lines are the two threshold levels. Both events clear **both**
thresholds — which is precisely why the second event is the interesting one.
"""),
code("""
t_h, tau_h = _s1_stress(HOURS)
u_h = stress_to_speed(tau_h)

fig, (a1, a2) = plt.subplots(2, 1, figsize=(9.5, 5.2), sharex=True,
                             gridspec_kw={'height_ratios': [2, 1]})
a1.fill_between(t_h, 0, tau_h, color=CLR[0], alpha=.30)
a1.plot(t_h, tau_h, color=CLR[0], lw=2)
for lvl, lab, c in [(0.05, r'$\\tau_{ce}$ mud = 0.05 Pa', CLR[2]),
                    (0.10, r'$\\tau_{ce}$ sand = 0.10 Pa', CLR[3])]:
    a1.axhline(lvl, color=c, lw=1.4, ls='--')
    a1.text(HOURS * .99, lvl * 1.10, lab, ha='right', fontsize=8.5, color=c)
a1.set_ylabel(r'$\\tau_b$  [Pa]'); a1.set_title('Prescribed bottom-stress history')
a1.annotate('event 1\\npeak 1.00 Pa', xy=(24, 1.0), xytext=(0, -34),
            textcoords='offset points', ha='center', fontsize=9,
            arrowprops=dict(arrowstyle='->', color='0.3'))
a1.annotate('event 2\\npeak 0.45 Pa', xy=(90, .45), xytext=(0, 30),
            textcoords='offset points', ha='center', fontsize=9,
            arrowprops=dict(arrowstyle='->', color='0.3'))
a2.plot(t_h, u_h, color='0.35', lw=1.8)
a2.set_ylabel('$|u|$ [m/s]'); a2.set_xlabel('time [h]')
a2.set_title(r'... imposed through the current, $\\tau_b=\\rho c_d |u|^2$', fontsize=9.5)
plt.tight_layout(); plt.show()
"""),
md("""
## Run

2000 elements, 500 per class, all starting **on the bed** (`z='seafloor'`,
`settled=1`). Five simulated days; this takes a few minutes.
"""),
code("""
ds = column_dataset(h=H, u=u_h, u_star=0.02, hours=HOURS, dt_hours=1.0, nz=81)
o, reader = column_model(
    ds,
    vertical_mixing__settling_model='prescribed',
    vertical_mixing__tau_crit_mode='constant',
    vertical_mixing__resuspension_height_mode='turbulent',
    vertical_mixing__bbl_scheme='legacy',
    vertical_mixing__resuspension_seed_layer=5.0)

np.random.seed(3)
n = N_PER_CLASS * 4
o.seed_elements(lon=8.0, lat=64.0, number=n, time=reader.start_time, z='seafloor',
                terminal_velocity=np.repeat(-CLASS_WS, N_PER_CLASS),
                tau_crit=np.repeat(CLASS_TAU_CE, N_PER_CLASS),
                grain_diameter=np.repeat(CLASS_D, N_PER_CLASS),
                rho_s=np.full(n, 2650.0), use_stokes=np.zeros(n),
                origin_marker=np.repeat(np.arange(4), N_PER_CLASS),
                settled=np.ones(n), moving=np.zeros(n))

out = os.path.join(tempfile.mkdtemp(), 's1.nc')
o.run(time_step=300, time_step_output=1800,
      duration=datetime.timedelta(hours=HOURS), outfile=out,
      export_variables=['z', 'settled', 'moving', 'origin_marker',
                        'times_resuspended'], stop_on_error=True)

with xr.open_dataset(out) as d:
    settled = d['settled'].values
    zz = d['z'].values
    marker = d['origin_marker'].values[:, 0].astype(int)
    resusp = d['times_resuspended'].values
hours = np.arange(settled.shape[1]) * 0.5
susp = np.array([[np.nanmean(settled[marker == k, i] == 0)
                  for i in range(settled.shape[1])] for k in range(4)])
print(f"lost from the domain: {int(np.sum(~np.isfinite(settled[:, -1])))} of {n}")
"""),
md("""
## Result 1 — suspended load and its termination

This is the central figure. Shaded bands are the two stress events.
"""),
code("""
fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 6.4), sharex=True,
                             gridspec_kw={'height_ratios': [1, 2.4]})
for a in (a1, a2):
    for t0, t1 in [(6, 42), (78, 102)]:
        a.axvspan(t0, t1, color='0.88', zorder=0)
a1.plot(t_h, tau_h, color='0.3', lw=1.8); a1.set_ylabel(r'$\\tau_b$ [Pa]')
a1.set_title('Suspended fraction by class')

for k, (nm, c) in enumerate(zip(NAMES, CLR)):
    a2.plot(hours, susp[k], color=c, lw=2.2, label=nm)
a2.set_xlabel('time [h]'); a2.set_ylabel('fraction in suspension')
a2.set_ylim(-.03, 1.05); a2.legend(loc='center right', fontsize=9)
a2.annotate('coarse sand deposits\\nfirst as stress falls',
            xy=(43, .05), xytext=(24, 40), textcoords='offset points',
            fontsize=9, arrowprops=dict(arrowstyle='->', color=CLR[3]))
a2.annotate('clay still suspended\\nafter 5 days', xy=(118, susp[0][-1]),
            xytext=(-135, 18), textcoords='offset points', fontsize=9,
            arrowprops=dict(arrowstyle='->', color=CLR[0]),
            bbox=dict(boxstyle='round,pad=0.35', fc='#eef5fa', ec=CLR[0]))
a2.annotate('event 2 lifts nearly as much\\nsand as event 1 — it should not',
            xy=(90, susp[3][np.argmin(abs(hours - 90))]), xytext=(-40, -78),
            textcoords='offset points', fontsize=9, ha='center',
            arrowprops=dict(arrowstyle='->', color=CLR[3]),
            bbox=dict(boxstyle='round,pad=0.35', fc='#fff3f0', ec=CLR[3]))
plt.tight_layout(); plt.show()
"""),
md("""
### Deposition order

The order in which classes clear from suspension after event 1 ends is the
Lagrangian expression of the fining-upward storm layer.
"""),
code("""
clear_h = []
for k in range(4):
    after = np.where((hours >= 42) & (susp[k] < 0.10))[0]
    clear_h.append(hours[after[0]] if after.size else np.inf)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.2))
lab = [f'{h:.1f} h' if np.isfinite(h) else 'never' for h in clear_h]
vals = [h if np.isfinite(h) else HOURS for h in clear_h]
a1.barh(NAMES, vals, color=[CLR[2] if np.isfinite(h) else CLR[0] for h in clear_h], alpha=.9)
a1.axvline(42, color='k', lw=1.4, ls='--')
a1.text(42.8, 3.35, 'event 1 ends', fontsize=8.5)
for i, (v, t) in enumerate(zip(vals, lab)):
    a1.text(v + 1.5, i, t, va='center', fontsize=9.5, fontweight='bold')
a1.set_xlabel('hour at which <10 % remains suspended')
a1.set_title('Deposition order follows settling velocity')

a2.loglog(CLASS_WS * 1e3, [max(h - 42, .2) for h in vals], 'o-',
          color=CLR[0], ms=8, lw=1.8)
for w, v, nm in zip(CLASS_WS * 1e3, vals, NAMES):
    a2.annotate(nm.split()[0], (w, max(v - 42, .2)), textcoords='offset points',
                xytext=(6, 6), fontsize=8.5)
a2.set_xlabel('settling velocity [mm/s]'); a2.set_ylabel('hours to clear after event 1')
a2.set_title('Faster settlers clear sooner')
plt.tight_layout(); plt.show()
for nm, t in zip(NAMES, lab):
    print(f'{nm:>14s}: {t}')
"""),
md("""
## Result 2 — where it fails

Two failures, sharing one cause: our resuspension is an **all-or-nothing
threshold crossing**, whereas ROMS erodes with a *flux*,
$E = E_0(1-\\phi)(\\tau_b/\\tau_{ce}-1)$, which scales with the excess stress.
"""),
code("""
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.2))

# --- F2: event-2 intensity is not suppressed
e1 = np.array([susp[k][(hours >= 6) & (hours <= 42)].max() for k in range(4)])
e2 = np.array([susp[k][(hours >= 78) & (hours <= 102)].max() for k in range(4)])
x = np.arange(4)
a1.bar(x - .19, e1, .38, color='0.35', label='event 1 (1.00 Pa)')
a1.bar(x + .19, e2, .38, color=CLR[3], alpha=.9, label='event 2 (0.45 Pa)')
a1.axhline(e1[3] * 0.25, color=CLR[2], lw=1.8, ls='--')
a1.text(3.42, e1[3] * 0.25 + .03, 'what Sherwood\\nimplies for sand',
        ha='right', fontsize=8.5, color=CLR[2])
a1.set_xticks(x); a1.set_xticklabels([n.split()[0] + '\\n' + n.split()[1] for n in NAMES],
                                     fontsize=8.5)
a1.set_ylabel('peak suspended fraction'); a1.set_ylim(0, 1.15)
a1.set_title('F2: a weaker event lifts nearly as much'); a1.legend(fontsize=8.5)
a1.text(3, e2[3] + .05, f'{e2[3]/e1[3]:.2f}×', ha='center', fontsize=10,
        fontweight='bold', color=CLR[3])

# --- F3: lift-offs per particle
per = []
for k in range(4):
    m = marker == k
    last = np.nanmax(np.where(np.isfinite(resusp[m]), resusp[m], np.nan), axis=1)
    per.append(np.nanmean(last))
a2.bar(NAMES, per, color=[CLR[2] if p < 5 else CLR[3] for p in per], alpha=.9)
a2.axhline(4, color=CLR[2], lw=1.8, ls='--')
a2.set_yscale('log')
a2.text(3.42, 4.8, 'a few per event\\nwould be physical', ha='right',
        fontsize=8.5, color=CLR[2])
a2.set_ylabel('mean lift-offs per particle  (log)')
a2.set_title('F3: settled grains re-lift every time step')
for i, p in enumerate(per):
    a2.text(i, p * 1.15, f'{p:.0f}', ha='center', fontsize=10, fontweight='bold')
plt.setp(a2.get_xticklabels(), fontsize=8.5)
plt.tight_layout(); plt.show()

print('peak suspended fraction   event 1 / event 2 / ratio')
for nm, a, b in zip(NAMES, e1, e2):
    print(f'  {nm:>14s}  {a:.3f}  {b:.3f}   {b/a:.2f}')
print('\\nmean lift-offs per particle:', ', '.join(f'{n.split()[1]} {p:.0f}' for n, p in zip(NAMES, per)))
"""),
md("""
## Verdict

**Four passes.**

| | Sherwood | model |
|---|---|---|
| Fines dominate the suspended load | yes | peak fraction falls monotonically with grain size |
| Coarse deposits first | yes | 140 µm → 62.5 → 30 µm, strictly ordered |
| 4 µm still suspended at day 5 | *"mostly"* | ~64 % |
| Mass conserved in the column | — | 0 elements lost |

**Two failures, one cause.**

- **F2 — event-intensity scaling.** The second event, at 0.45 Pa against a 0.10 Pa
  threshold, resuspends ~0.92× as much sand as the 1.00 Pa first event. Sherwood
  gets *"minimal amounts"* because his erosion is a flux proportional to the
  excess stress. A threshold crossing has no notion of *how far* it was exceeded.

- **F3 — hopping.** A settled element above threshold is lifted **again on the
  next time step**, and the coarse classes — which fall back quickly — cycle
  continuously: **~299 lift-offs per sand particle** over two events, where a
  handful would be physical.

  This diagnostic was itself broken until this campaign: `times_resuspended` was
  declared `uint8` and wrapped silently at 255, which is why the same quantity
  came out as 44 in one analysis and 254 in another. Neither was real. It is now
  `uint32` — and note that the counter is exported by the DDT production runs, so
  resuspension counts in existing output are wrapped wherever they exceeded 255.

  **F3 matters beyond this notebook.** Each hop displaces a particle downstream by
  a step of advection that has nothing to do with the physics. This is the
  mechanism behind the parent DDT study's ~45 % out-of-domain loss of its fastest
  settling classes, which currently blocks the class-weight inversion — isolated
  here in a case with a known answer.

**Consequence.** Probabilistic erosion (an Ariathurai–Partheniades pickup flux) is
**required**, not an optional refinement. M1 hands that to **M2** as its first
item. A third gap surfaced along the way: Sherwood also has a critical shear
stress for *deposition* ($\\tau_d$, Krone 1962) that we lack entirely — we deposit
unconditionally.
"""),
    ])


# =========================================================================
@builder("00_overview.ipynb")
def _overview():
    return nb([
md("""
# SPPMDRIFT — M1 resuspension validation

This is the index for the M1 validation campaign. Each case has its own notebook;
this one summarises what was tested, what passed, and what did not.

**Why this exists.** The resuspension physics in this fork — a Styles & Glenn
(2000) wave–current bed stress, a dynamic critical shear stress with fractal floc
strength and consolidation, and a Rouse-based lift-off height — was already
built, and has been the production default since June 2026. What it had never
had was evidence that it is *right*. Before that physics is carried into
SPPMDRIFT, it needs benchmarks with known answers.

**The standing rule.** Physics is validated against **published idealized cases**,
never against the parent study's observed sediment data. An observed deposition
pattern confounds transport physics with source history, so it cannot isolate a
physics defect.

## The notebooks

| | Case | Tests | Verdict |
|---|---|---|---|
| [01](01_settling_column.ipynb) | Settling column | settling in isolation | **pass** |
| [02](02_well_mixed_condition.ipynb) | Well-mixed condition | the vertical random walk | **pass** |
| [03](03_rouse_profile.ipynb) | Rouse equilibrium profile | settling + mixing + lift-off height, jointly | **pass**, with a documented bias |
| [04](04_S3_mixed_bed_closure.ipynb) | Sherwood Fig. 2 — mixed-bed threshold | the erosion threshold of a sand–mud bed | **fail (F1)** |
| [05](05_S1_double_resuspension.ipynb) | Sherwood §3.2.1 — double resuspension | erosion and deposition through two storm events | **4 pass, 2 fail (F2, F3)** |

**Tier A** (01–03) are analytic: closed-form answers, no external data, seconds to
run. **Tier B** (04–05) are adapted from the ROMS/COAWST sediment model papers:

- Warner, J. C. *et al.* (2008), *Computers & Geosciences* **34**, 1284–1306.
- Sherwood, C. R. *et al.* (2018), *Geosci. Model Dev.* **11**, 1849–1871.

## The verdict in one line

**The settling, mixing and deposition chain is sound. The erosion side is not.**
"""),
header("M1 validation - overview"),
md("""
## Scoreboard

Every case, its metric and its outcome. Green is a pass; red is a defect that is
pinned by a test rather than tuned away.
"""),
code("""
rows = [
    # tier, case, metric, reference, model, ok
    ('A', 'Settling column',        'settled fraction @2 h',      0.720, 0.719, True),
    ('A', 'Well-mixed condition',   'max interior deviation',     0.00,  0.09,  True),
    ('A', 'Rouse profile P=0.25',   'fitted exponent',            0.25,  0.223, True),
    ('A', 'Rouse profile P=0.50',   'fitted exponent',            0.50,  0.404, True),
    ('B', 'S3 low-mud limit',       'tau_ce / Shields',           1.00,  1.000, True),
    ('B', 'S3 cohesive limit 4 um', 'ours / Sherwood Eq. 6',      1.00,  0.6,   False),
    ('B', 'S3 cohesive limit 140um','ours / Sherwood Eq. 6',      1.00,  14.1,  False),
    ('B', 'S1 suspended ordering',  'peak fraction, 140 um',      0.00,  0.66,  True),
    ('B', 'S1 fines @ day 5',       'suspended fraction, 4 um',   0.50,  0.64,  True),
    ('B', 'S1 event-2 suppression', 'peak2 / peak1, 140 um',      0.25,  0.92,  False),
    ('B', 'S1 lift-offs',           'per particle, 140 um',       4.00,  299.0, False),
    ('C', 'SG2000 Table 2',         'max abs error (15 entries)', 0.00,  0.000, True),
    ('C', 'Settling vs Maggi 2013', 'log10-RMSE, bb16',           0.00,  0.202, True),
]
labels = [f"{t} · {c}" for t, c, *_ in rows]
ok = [r[-1] for r in rows]

fig, ax = plt.subplots(figsize=(9.5, 5.6))
y = np.arange(len(rows))[::-1]
ax.barh(y, [1] * len(rows),
        color=['#e8f3ec' if o else '#fdeceb' for o in ok], height=.82)
for yi, (tier, case, metric, ref, mod, o) in zip(y, rows):
    ax.text(0.015, yi, f"{tier} · {case}", va='center', fontsize=9.5,
            fontweight='bold' if not o else 'normal')
    ax.text(0.50, yi, metric, va='center', fontsize=8.5, color='0.35')
    ax.text(0.885, yi, f"{mod:g}", va='center', ha='right', fontsize=9.5,
            fontweight='bold', color='#1d6b3a' if o else '#a8302c')
    ax.text(0.985, yi, 'PASS' if o else 'FAIL', va='center', ha='right',
            fontsize=9, fontweight='bold', color='#1d6b3a' if o else '#a8302c')
ax.set_xlim(0, 1); ax.set_ylim(-.6, len(rows) - .4)
ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
for s in ax.spines.values():
    s.set_visible(False)
ax.set_title('M1 validation scoreboard  —  %d pass, %d fail'
             % (sum(ok), len(ok) - sum(ok)), fontsize=12, pad=12)
plt.tight_layout(); plt.show()
"""),
md("""
## The four findings

All four are consequences of how erosion is formulated, and the last three share
a single root cause.

### F1 — the mixed-bed threshold is not Sherwood's

Our `mixed` mode blends the particle's **own** fractal-floc strength, which scales
with its own diameter, instead of the **bulk bed** stress, and it omits the
`max(·, τ_c)` floor. Errors run from 0.6× (clay, too easy to erode) to 14×
(coarse sand, far too hard). → notebook [04](04_S3_mixed_bed_closure.ipynb)

### F2 — a threshold cannot scale with event intensity

ROMS erodes with a flux proportional to the excess stress. We lift *everything*
above threshold, however small the excess, so a weaker second storm resuspends
0.92× as much sand as a storm twice its size. → notebook
[05](05_S1_double_resuspension.ipynb)

### F3 — settled grains re-lift on every time step

**299** lift-offs per sand particle across two storms, where a handful would be
physical. (This count was itself corrupted by a `uint8` counter that wrapped at
255, found and fixed during the campaign — see F6 in the report.) **This is the mechanism behind the parent DDT study's ~45 % loss of its
fastest settling classes out of the domain**, which currently blocks the
class-weight inversion. → notebook [05](05_S1_double_resuspension.ipynb)

### F4 — no critical shear stress for deposition

Sherwood implements Krone's $\\tau_d$: deposition is suppressed while the bed
stress is high. We deposit unconditionally whenever an element reaches the bed.
Not anticipated when M1 was planned; surfaced by reading the reference.

## What follows

F2 and F3 share a cause, so **probabilistic erosion is required, not optional**.
M1 hands that to M2 as its first item, together with F1 and F4. The passes in
Tier A mean the settling and mixing machinery underneath can be trusted while
that work happens.

Full write-up: [`../REPORT.md`](../REPORT.md). Module plan:
[`../../M1_resuspension_validation.md`](../../M1_resuspension_validation.md).
"""),
    ])


if __name__ == "__main__":
    main()
