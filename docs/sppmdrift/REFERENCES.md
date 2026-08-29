# SPPMDRIFT — seed bibliography

**Read this first.** The `Status` column is about *provenance*, not quality:

| Status | Meaning |
|---|---|
| `read` | The paper is on disk and has been read for this project. |
| `on disk` | The PDF/HTML is present; the citation reflects its filename or the code that cites it. Details still worth confirming. |
| `unvetted` | **Proposed from memory.** The work is real and relevant, but author list, year, volume and page numbers have **not** been checked against the paper. Treat every such line as a lead, not a citation. |

Nothing here has been through a module's §3 gate except the two ROMS/COAWST
papers used by M1. Pruning and extending these lists *is* the gate.

Local copies live outside this repository, in `ddt_dump/resuspension_mechanics/`
and `ddt_dump/bblm_sg2000/`.

---

## Implemented physics — the papers the current code is built on

| Ref | Used by | Status |
|---|---|---|
| Styles, R. & Glenn, S. M. (2000). Modeling stratified wave and current bottom boundary layers on the continental shelf. *J. Geophys. Res.* 105(C10), 24119–24139. | `bblm_sg2000.py` | on disk |
| Styles, Glenn & Brown (2017). ERDC/CHL TR-17-11 (arbitrary roughness). | `bblm_sg2000.py` | unvetted |
| Wiberg, P. L. & Smith, J. D. (1987). Calculations of the critical shear stress for motion of uniform and heterogeneous sediments. | `floc_strength.py` | on disk |
| Soulsby, R. L. & Whitehouse, R. J. S. (1997). Threshold of sediment motion in coastal environments. | `tau_crit_shields` | unvetted |
| Kranenburg, C. (1994). The fractal structure of cohesive sediment aggregates. | `tau_crit_floc` | unvetted |
| Sanford, L. P. (2008). Modeling a dynamically varying mixed sediment bed with erosion, deposition, bioturbation, consolidation, and armoring. | `consolidate_phi` | unvetted |
| Maggi, F. (2013). The settling velocity of mineral, biomineral, and biological particles and aggregates in water. *J. Geophys. Res. Oceans* 118. | `settling.py`, test fixture | on disk (data) |
| Bagheri, G. & Bonadonna, C. (2016). On the drag of freely falling non-spherical particles. | `settling.py` BB16 | unvetted |
| Dietrich, W. E. (1982). Settling velocity of natural particles. | `settling.py` | unvetted |
| Jacobs, W. et al. (2011). Erosion threshold of sand–mud mixtures. | `mixed` mode ramp | on disk |
| Yao, P. et al. (2022). Erosion behavior of sand–silt mixtures. | `mixed` mode ramp | on disk |
| Rouse, H. (1937). Modern conceptions of the mechanics of fluid turbulence. | `rouse.py` | unvetted |

## M1 — resuspension validation

**Gate cleared.** The two ROMS/COAWST papers below are `read` and vetted for this
purpose; the analytic references are not.

| Ref | Role | Status |
|---|---|---|
| **Warner, J. C., Sherwood, C. R., Signell, R. P., Harris, C. K. & Arango, H. G. (2008). Development of a three-dimensional, regional, coupled wave, current, and sediment-transport model. *Computers & Geosciences* 34, 1284–1306.** | W1 open-channel, W2 migrating trench | **read** |
| **Sherwood, C. R., Aretxabaleta, A. L., Harris, C. K., Rinehimer, J. P., Verney, R. & Ferré, B. (2018). Cohesive and mixed sediment in ROMS (v3.6) implemented in COAWST. *Geosci. Model Dev.* 11, 1849–1871.** | S1, S2, S3 | **read** |
| van Rijn, L. C. (1987/1993). Trench flume experiment. | W2 source data | unvetted |
| Thomson, D. J. (1987). Criteria for the selection of stochastic models of particle trajectories in turbulent flows. | well-mixed condition | unvetted |
| Visser, A. W. (1997). Using random walk models to simulate the vertical distribution of particles in a turbulent water column. | well-mixed condition | unvetted |
| Ross, O. N. & Sharples, J. (2004). Recipe for 1-D Lagrangian particle tracking models in space-varying diffusivity. | well-mixed condition | unvetted |
| Grant, W. D. & Madsen, O. S. (1979). Combined wave and current interaction with a rough bottom. | BBL pure-wave limit | unvetted |
| Krone, R. B. (1962). Flume studies of the transport of sediment in estuarial shoaling processes. | deposition benchmark | unvetted |
| Partheniades, E. (1965). Erosion and deposition of cohesive soils. | erosion benchmark | unvetted |
| van Rijn, L. C. (1984). Sediment transport, part II: suspended load transport. | reference concentration | unvetted |
| Sanford, J. P. & Maa, J. P.-Y. (2001). A unified erosion formulation for fine sediments. *Marine Geology*. | erosion benchmark | on disk (HTML) |
| Ariathurai, R. & Arulanandan, K. (1978). Erosion rates of cohesive soils. | pickup flux | on disk |
| Henry, C. et al. (2023). Particle resuspension: challenges and perspectives. | review | on disk |

## M2 — burial and thresholds  *(gate open; all unvetted unless noted)*

Sanford (2008) — the central reference; Sanford & Maa (2001) *(on disk)*;
Harris, C. K. & Wiberg, P. L. (2001), suspended transport and bed reworking on
shelves; Le Hir, P., Cayocca, F. & Waeles, B. (2011), erosion of cohesive and
mixed sediment; Winterwerp, J. C. & van Kesteren, W. G. M. (2004), *Introduction
to the Physics of Cohesive Sediment in the Marine Environment* (book);
Winterwerp et al. (2012) *(on disk)*; Mehta, A. J. (2013), *An Introduction to
Hydraulics of Fine Sediment Transport*.

Bioturbation: Boudreau, B. P. (1998), Mean mixed depth of sediments: the
wherefores and the whys — source of the ~9.8 cm global mean mixed-layer depth;
Boudreau (1986), diffusive models of bioturbation; Meysman, F. J. R.,
Middelburg, J. J. & Heip, C. H. R. (2006), Bioturbation: a fresh look at Darwin's
last idea; Teal, L. R. et al. (2008), global patterns of bioturbation intensity
and mixed depth; Wheatcroft, R. A. et al. (1990), tracer burial.

Armoring: Wiberg, P. L., Drake, D. E. & Cacchione, D. A. (1994), sediment
resuspension and bed armoring during high bottom stress events *(cited by
Sherwood et al. 2018)*.

Pollutant-specific: Eganhouse, R. P. & Pontolillo, J. (2008), DDE in sediments of
the Palos Verdes shelf — in-situ transformation rates and geochemical fate.

## M3 — settling  *(gate open; all unvetted)*

Ferguson, R. I. & Church, M. (2004), a simple universal equation for grain
settling velocity; Dioguardi, F., Mele, D. & Dellino, P. (2018), drag for
irregularly shaped particles over a wide Reynolds range; Corey, A. T. (1949),
shape factor; Khelifa, A. & Hill, P. S. (2006), models for effective density and
settling velocity of flocs — **the named variable-`Df` gap**;
Winterwerp, J. C. (1998, 2002); Strom, K. & Keyvani, A. (2011), an explicit
full-range settling velocity equation for mud flocs; Manning, A. J. & Dyer, K. R.
(2007); Waldschläger, K. & Schüttrumpf, H. (2019), settling and rise of
microplastics; Khatmullina, L. & Isachenko, I. (2017); Kooi, M. et al. (2017),
biofouling and the vertical transport of microplastics.

## M4 — chemistry and the dissolved fraction  *(gate open; all unvetted)*

The model: Aghito, M., Calgaro, L., Dagestad, K.-F., Ferrarin, C., Marcomini, A.,
Breivik, Ø. & Hole, L. R. (2023), ChemicalDrift 1.0, *Geosci. Model Dev.* 16;
Dagestad, K.-F., Röhrs, J., Breivik, Ø. & Ådlandsvik, B. (2018), OpenDrift v1.0,
*Geosci. Model Dev.* 11.

Partitioning: Karickhoff, S. W., Brown, D. S. & Scott, T. A. (1979), sorption of
hydrophobic pollutants on natural sediments; Karickhoff & Morris (1985) *(both
already used inside `chemicaldrift.py`)*; Park & Clough (2014) *(the KOC
regressions in the code)*; Schwarzenbach, R. P., Gschwend, P. M. & Imboden, D. M.,
*Environmental Organic Chemistry* (textbook); Burkhard, L. P. (2000), DOC
partitioning; Accardi-Dey, A. & Gschwend, P. M. (2002), natural organic matter and
black carbon as sorbents in sediments.

Sorption kinetics: Pignatello, J. J. & Xing, B. (1996), slow sorption;
Cornelissen, G. et al. (1997, 1998, 2005), rapidly and slowly desorbing fractions.

Sediment–water exchange: Boudreau, B. P. (1997), *Diagenetic Models and Their
Implementation*; Thibodeaux, L. J., *Environmental Chemodynamics*.

Mass-transfer particle tracking — the formulation anchor for a per-element mass
vector: Benson, D. A. & Bolster, D. (2016), arbitrarily complex chemical reactions
on particles; Schmidt, M. J., Pankavich, S. & Benson, D. A. (2017), a kernel-based
Lagrangian method for imperfectly mixed reactions; Engdahl, N. B., Benson, D. A. &
Bolster, D. (2017).

The compound: Eganhouse & Pontolillo (2008), as above.

## M5 — dynamic particle size  *(gate open; all unvetted)*

Winterwerp, J. C. (1998), a simple model for turbulence-induced flocculation of
cohesive sediment — the single-class kinetic equation; Winterwerp (2002);
Khelifa & Hill (2006); Keyvani, A. & Strom, K. (2014), floc growth and breakup
under cycled shear *(a Sherwood et al. 2018 §3.1.1 validation case)*;
Verney, R., Lafite, R., Brun-Cottan, J.-C. & Le Hir, P. (2011), floc population
over a tidal cycle — FLOCMOD *(the other §3.1.1 case)*; Lee, B. J., Toorman, E.,
Molz, F. J. & Wang, J. (2011), a two-class population balance equation;
Maerz, J. et al. (2011); Smoluchowski, M. (1917); Jackson, G. A. (1990), marine
algal floc formation by physical coagulation; Burd, A. B. & Jackson, G. A. (2009),
Particle aggregation, *Annual Review of Marine Science* 1 *(review)*.

## Numerics (cross-cutting)  *(unvetted)*

Nordam, T. & Duran, R. (2020), numerical integrators for Lagrangian oceanography;
Gräwe, U. (2011), high-order particle-tracking schemes in a water column model.

## Explicitly excluded from physics validation

**Wu, M. S. C., Schmidt, J. T., Kittner, H. E., Earth 182B Group & Valentine,
D. L. (2025). Geospatial and Temporal Inventory for Industrial DDT Waste Disposal
to a Deep Coastal Ocean Environment. *Environ. Sci. Technol.* 59, 20578–20587.
doi:10.1021/acs.est.5c03851** — *on disk*.

This is the **science** target of the parent DDT study (a deposition pattern), not
a physics benchmark. An observed deposition pattern confounds transport physics
with source history and cannot isolate a physics defect. Do not use it to validate
physics. See `CLAUDE.md`.
