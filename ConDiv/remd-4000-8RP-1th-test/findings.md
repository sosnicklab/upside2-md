# Findings

- 2026-02-12: `run.output` contains mixed historical runs; one logged crash is `ModuleNotFoundError: No module named 'numpy'` from `ConDiv.py`, consistent with launching via system `python3` (`/opt/homebrew/opt/python@3.14/bin/python3.14`) instead of the project venv interpreter.
- 2026-02-12: `ConDiv.py` had duplicated `__main__` dispatch (`worker`/`restart`/`initialize`) not present in `ConDiv_original.py`, causing mode handlers to execute twice.
- 2026-02-12: Existing checkpoint worker script path is `test_00/ConDiv.py`; restart runs use that copy, so fixes to top-level `ConDiv.py` must be mirrored there to affect worker behavior.
- 2026-02-12: RCSB endpoints selected for automation:
  - Search API: `https://search.rcsb.org/rcsbsearch/v2/query`
  - Entry metadata API: `https://data.rcsb.org/rest/v1/core/entry/{pdb_id}`
  - Coordinate download: `https://files.rcsb.org/download/{pdb_id}.pdb`
- 2026-02-12: Current training-set-like filters inferred from project data and code:
  - sequence-length range in existing `upside_input` is 50 to 151 residues
  - chain-break exclusion criterion from `ConDiv.py`: `max_sep < 2.0`
- 2026-02-12: In `ConDiv_original.py` training-set acceptance logic, the actual structural filter is applied after loading converted coordinates: keep only proteins with `max_sep < 2.0` (no chain breaks). No additional explicit biochemical/metadata filters are present there.
