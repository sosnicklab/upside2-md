2026-03-09
- Initialized task tracking files for the `ConDiv_symlay` workflow investigation.
- Traced `ConDiv_symlay/run_init.sh`, `ConDiv_mem.py`, `run_upside.py`, `upside_config.py`, `src/membrane_potential.cpp`, and `symlay_utils.py` to identify the active membrane potential form and training constraints.
- Verified with local inspection that no initialized `ConDiv_symlay/test_dimer3` run is present in the repository; only source `param*/membrane.h5` files and a test layer manifest are checked in.
- Ran local read-only validation against `ConDiv_symlay/param0/membrane.h5` and `ConDiv_symlay/test_symlay_manifest/layer_manifest.json`; confirmed the checked-in source membrane is not yet in the symmetric constrained subspace and inferred the required projected supports from the manifest.
- Traced the reported HDF5/open-config failure through the `ConDiv_symlay` worker path and identified a nested-`srun` launch bug: the Python worker was already running inside a Slurm step, but it launched the non-MPI `upside` binary through a second `srun --ntasks=n_replica`.
- Patched `ConDiv_symlay/ConDiv_mem.py` so the worker launches `upside` locally inside the existing worker step, and so the temporary template config is moved into `run.0.h5` instead of being redundantly regenerated and retained as `base.h5`.
- Updated `ConDiv_symlay/README.md` to document the corrected Slurm launch model.
- Verified locally that `ConDiv_symlay/ConDiv_mem.py` compiles, that the revised config-generation path produces valid `run.0.h5` and `run.1.h5`, and that the temporary `base.h5` is no longer left behind in the smoke directory.
