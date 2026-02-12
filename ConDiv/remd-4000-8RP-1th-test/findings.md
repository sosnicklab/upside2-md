# Findings

- 2026-02-12: `run.output` contains mixed historical runs; one logged crash is `ModuleNotFoundError: No module named 'numpy'` from `ConDiv.py`, consistent with launching via system `python3` (`/opt/homebrew/opt/python@3.14/bin/python3.14`) instead of the project venv interpreter.
- 2026-02-12: `ConDiv.py` had duplicated `__main__` dispatch (`worker`/`restart`/`initialize`) not present in `ConDiv_original.py`, causing mode handlers to execute twice.
- 2026-02-12: Existing checkpoint worker script path is `test_00/ConDiv.py`; restart runs use that copy, so fixes to top-level `ConDiv.py` must be mirrored there to affect worker behavior.
