# Scope boundary — read before working in this directory

This directory is the **SPPMDRIFT model-development project**. It is a separate
offspring of the DDT dispersal study in `~/projects/science/ddt_dump`, and the
two must not be conflated.

## If you are working on the DDT study (`ddt_dump`)

Do **not** pull this directory in. Nothing here is needed to run, analyse or
write up the DDT simulations. The DDT study's own design documents
(`WORKFLOW_PLAN.md`, `BBL_SG2000_PLAN.md`, `RESUSPENSION_MECHANICS_PLAN.md`,
`SETTLING_IMPLEMENTATION.md`, `RECAP.md`, `TODO.md`) remain the authority on what
was built there and why.

## If you are working on SPPMDRIFT

- Treat `ddt_dump` as a **read-only** source of data, forcing, papers and prior
  results. Do not modify it, and do not re-derive its science conclusions.
- `SedimentDrift` must stay **behaviour-identical**. DDT production brackets
  (11 velocity classes x 162k particles per scenario) are reproducible only if it
  does. New physics goes behind a config default that preserves current
  behaviour, or into `SPPMDrift`.
- Respect the **separation-of-concerns contract** in `ARCHITECTURE.md`. It is the
  rule that keeps bed, particle, chemistry and transport concerns from tangling.
- **Do not start a module's design or code before its §3 literature gate has been
  cleared jointly with S. Soares.** The candidate citations in the module docs
  are unvetted proposals.

## Things that are settled — do not relitigate

- Plans live here, in the fork, version-controlled beside the code.
- Architecture is *new class + shared physics package*, not "evolve
  sedimentdrift.py in place" and not "subclass SedimentDrift".
- Wu et al. (2025) sediment DDX is the DDT study's **science** target. It is
  **not** a physics-validation benchmark and must not be used as one.
- Physics validation uses published idealized / semi-idealized cases —
  Warner et al. (2008) and Sherwood et al. (2018) for the ROMS/COAWST comparison.
- The dissolved-phase spatial representation (M4) is **deliberately deferred** to
  that module's design session. Three candidates are recorded; none is chosen.
