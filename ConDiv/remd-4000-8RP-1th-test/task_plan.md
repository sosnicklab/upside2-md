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
- Add transient-network hardening to RCSB metadata/download calls so chunked transfer interruptions do not crash batch runs.

## Active Debug Task (2026-02-14)

### Project Goal
Debug `prepare_rcsb_upside_input.py` crash caused by `http.client.IncompleteRead` during PDB download and keep batch processing resilient.

### Architecture & Key Decisions
- Keep the existing pipeline shape (fetch metadata -> download PDB -> convert/filter); do not alter dataset filtering semantics.
- Add bounded retries with short backoff around network operations and treat incomplete/chunked read interruptions as retryable transport failures.
- Preserve per-entry skip behavior: when retries are exhausted, mark the entry as skipped with a specific reason and continue.

### Execution Phases
- [x] Phase A: Confirm failing call path and identify unhandled exceptions.
- [x] Phase B: Implement retry/backoff for network reads (entry fetch + PDB download).
- [x] Phase C: Validate script syntax and log the fix.

### Known Errors / Blockers
- Live end-to-end validation is still limited by this environment's network restrictions.

## Active Debug Task (2026-02-16)

### Project Goal
Debug repeated `TypeError: model must be an integer, 1 is invalid` during NMR conversion and determine why accepted NMR count is zero.

### Architecture & Key Decisions
- Keep dataset filter semantics unchanged; isolate fix to argument typing at the converter boundary.
- Enforce integer typing for model selectors in both caller (`prepare_rcsb_upside_input.py`) and callee (`py/PDB_to_initial_structure.py`) to prevent string propagation.
- Use existing `manifest.csv` reason codes to identify where NMR entries were filtered out in the current run.

### Execution Phases
- [x] Phase A: Trace failure path for `--model` from NMR branch to ProDy parser.
- [x] Phase B: Apply minimal typing fixes for model arguments.
- [x] Phase C: Validate on local previously failing NMR samples and summarize filter-stage breakdown from manifest.

### Known Errors / Blockers
- `rcsb_dataset_all/manifest.csv` and list outputs were generated before this fix; a fresh full run is required to produce updated accepted-NMR counts.
