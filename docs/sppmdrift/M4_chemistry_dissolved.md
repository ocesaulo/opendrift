# M4 — Chemistry and the dissolved fraction

## 1. Goal & scope

Port the bed and carrier physics of `SedimentDrift` into the `ChemicalDrift`
toolkit, so that a tracked packet of pollutant can move between **particulate and
dissolved** form — fractioning its mass rather than switching identity — while
the sediment bed and the water column are genuinely decoupled.

This is the module in which SPPMDRIFT becomes itself: the element stops being a
grain of sediment and becomes **the centre of mass of a pollutant**.

Out of scope: dynamic particle size (M5).

## 2. Current state

**Nothing exists.** `ChemicalDrift` is never mentioned in any DDT planning
document, and `opendrift/models/chemicaldrift.py` (3019 lines) is stock upstream,
untouched by this fork. The DDT project has no chemistry, sorption or dissolved
phase anywhere.

`ChemicalDrift` as it stands:

- **Phases as an integer.** Each element has a `specie` index into up to 12
  species (`LMM` dissolved, `Humic colloid` = DOC, `Particle reversible` = SPM,
  `Sediment reversible`, `Sediment slowly reversible` = the buried pool, ...) and a
  scalar `mass` in µg. Mass is never split between phases.
- **Phase change is a stochastic jump.** `update_partitioning` (`:1415`) computes
  `p = 1 - exp(-K dt)`, draws a random number per element, and moves the whole
  element to a new specie — via a Python `for` loop over every transforming
  element (`:1443`).
- **Chemistry is good.** `init_transfer_rates` (`:687`) builds an
  `nspecies x nspecies` rate matrix from `LogKOW`, with `KOC_DOM = 2.88 KOW^0.67`
  and `KOC_sed = 2.62 KOW^0.82`, dissociation corrections from pKa and pH, and
  temperature/salinity corrections applied per step in `update_transfer_rates`.
- **Bed physics is far behind `SedimentDrift`.** `resuspension` (`:1607`) is a
  bare depth-averaged current-speed test against `resuspension_critvel` (default
  0.01 m/s) — no bed shear stress, no BBL, no `tau_crit`. `update_terminal_velocity`
  (`:1114`) is a hard-coded Sundby Stokes law on `diameter` and `density`.
  Burial is a single matrix entry, `[srev -> ssrev] = burial_rate / mixing_depth`,
  built **once at seed time** and never state-dependent. The code carries its own
  TODO that buried sediment ought to be a separate specie.

So the port is **bidirectional**, and neither direction is a copy.

## 3. Literature gate — OPEN

**Blocking.** Questions:

1. **Mass vector or specie index?** Is replacing the stochastic specie jump with a
   per-element mass vector across phases defensible, and what does the
   particle-tracking literature say about the resulting numerics? (Argued in
   `ARCHITECTURE.md`; the gate must confirm or refute it.)
2. **What velocity does a mixed-phase element move at?** The deferred decision —
   see §7. What is the criterion for when the choice matters?
3. **Equilibrium or kinetic sorption?** `ChemicalDrift` already uses a kinetic
   `k_ads`/`k_des`. For strongly hydrophobic legacy compounds, does slow
   desorption from a second, resistant domain dominate on decadal timescales?
4. **What carries the pollutant?** Organic carbon, or also black carbon/soot,
   which for DDT-class compounds can dominate sorption by an order of magnitude?
5. **Bed–water exchange.** What sets the flux out of a buried deposit — pore-water
   diffusion, bioirrigation, resuspension? This is what "decoupling sediment and
   water column" has to mean quantitatively.
6. **Compound properties.** What are defensible `LogKOW`, degradation half-lives
   and Henry constants for DDT, DDE and DDD? None are in
   `init_chemical_compound`'s table of 26 PAHs and 7 metals.

Candidate papers — **all unvetted**:

*The model to port into*
- Aghito et al. (2023), ChemicalDrift 1.0, *Geosci. Model Dev.*
- Dagestad et al. (2018), OpenDrift v1.0, *Geosci. Model Dev.*

*Partitioning*
- Karickhoff et al. (1979); Karickhoff & Morris (1985) *(already used in the code)*
- Schwarzenbach, Gschwend & Imboden, *Environmental Organic Chemistry* *(textbook)*
- Burkhard (2000), DOC partitioning
- Accardi-Dey & Gschwend (2002), black carbon as a sorbent

*Sorption kinetics / resistant desorption*
- Pignatello & Xing (1996); Cornelissen et al. (1997, 1998, 2005)

*Sediment–water exchange*
- Boudreau (1997), *Diagenetic Models and Their Implementation*
- Thibodeaux, *Environmental Chemodynamics*

*Mass-transfer particle tracking — the formulation anchor*
- Benson & Bolster (2016), arbitrarily complex chemical reactions on particles
- Schmidt, Pankavich & Benson (2017), a kernel-based Lagrangian method
- Engdahl, Benson & Bolster (2017)

*The compound*
- Eganhouse & Pontolillo (2008), DDE on the Palos Verdes shelf

## 4. Design

*Stub — pending §3.* Two elements of it are already argued and should be tested,
not assumed:

- the per-element **mass vector** replacing `specie` (`ARCHITECTURE.md`);
- **pore water as its own phase**, so dissolved pollutant can leave a buried
  deposit without the particle moving — the coupling point to M2's `burial_depth`
  (deeper implies a longer diffusion path).

## 5. Implementation steps

*Stub — pending §3.*

## 6. Validation

`tests/models/test_chemicaldrift.py::test_chemicaldrift_partitioning_organics` is
the existing template — it runs 300 days of 200 trajectories and asserts on the
resulting partition between species. Two properties are worth adding regardless
of the design chosen: **mass conservation** across all phases plus the degraded
and volatilized sinks, and **agreement with the analytic two-phase equilibrium**
when rates are held constant. If the mass-vector formulation is adopted, it should
reproduce the stochastic scheme's ensemble mean as particle count grows — that is
a direct, falsifiable check.

## 7. Open questions & decisions log

- **DEFERRED BY DECISION (2026-08-29): how the dissolved fraction is carried
  spatially.** Three candidates, none chosen:
  1. *Dual-position element* — a particulate position plus a dissolved offset
     advected with the water parcel; exchange falls off with separation; split
     into a real second element past a threshold. Closest to the
     "particle-nearby-dissolution-parcel pair" idea.
  2. *Paired elements* — one particulate and one dissolved element per packet,
     linked by an id. Reuses existing machinery; doubles particle count.
  3. *Single centre-of-mass element* — velocity is the mass-weighted mean of the
     settling and water velocities. Cheapest; cannot represent a plume separating
     into a settling limb and a dissolved limb.

  The criterion for when the choice matters: the separation accumulated over a
  phase-exchange timescale, `Delta ~ |w_s| tau_exchange`, against the plume scale.
  For a strongly-sorbing compound the dissolved fraction is small and (3) may
  suffice; for a weakly-sorbing one it will not. **This is measurable before
  committing, and should be measured.**
- Does the mass vector make `origin_marker`-style bookkeeping and the existing
  netCDF density-map writers inconsistent? `write_netcdf_chemical_density_map`
  assumes one specie per element.

## 8. Status

**Planned. Literature gate open — do not start §4.**
