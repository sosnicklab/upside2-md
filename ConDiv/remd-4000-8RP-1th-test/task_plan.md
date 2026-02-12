# Project Goal
Add a script that filters and downloads structures from RCSB, separates PDB/X-ray and NMR entries, excludes AlphaFold/ModelArchive-style computed structures, and prepares outputs in the same input format used under `upside_input/`.

# Architecture & Key Decisions
- Reuse existing project conventions for training inputs by inspecting current `upside_input/` and helper scripts, rather than inventing new formats.
- Implement a standalone preparation script that can:
  - query/filter RCSB entries,
  - split X-ray vs NMR lists,
  - download coordinate files,
  - run existing Upside conversion tooling to generate `.fasta`, `.initial.pkl`, and `.chi`.
- Keep filters configuration-driven via command-line arguments with safe defaults aligned to current dataset assumptions.

# Execution Phases
- [x] Phase 1: Inspect existing dataset format and filter logic in project scripts
- [x] Phase 2: Design RCSB filter/download/prep workflow and script interface
- [x] Phase 3: Implement script and output structure separation (xray/nmr)
- [x] Phase 4: Validate script usage and document run examples

# Known Errors / Blockers
- Full live end-to-end API/download run could not be completed in this environment due DNS/network restriction; script was syntax/CLI validated only.

# Revised Decisions
- Keep prior runtime fixes intact; extend repository with a new data-prep script for future dataset expansion.
