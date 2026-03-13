This example exercises the first ligand-aware Upside smoke workflow using the
4W52 T4 lysozyme L99A benzene complex.

What it does:
- prepares a protein-only Upside config with the existing tools
- patches ligand atoms and ligand interaction nodes directly into the `.up` file
- imports ligand and protein parameters from local CHARMM-format `.rtf`/`.prm` files
- adds example-local ligand-to-pocket distance restraints around the starting bound pose for smoke testing
- runs a short deterministic simulation
- validates that the restrained ligand-to-pocket geometry stays near the starting bound pose
- exports ligand-aware visualization files

Files:
- `0.prepare.py`: extracts chain A and benzene from `pdb/4W52.pdb`, builds the
  protein config, imports local CHARMM-format parameter files, then injects ligand
  nodes, starting-pose distance restraints, and pocket metadata
- `1.run.py`: runs a short smoke simulation
- `2.validate.py`: checks the trajectory and restrained ligand-to-pocket anchor geometry
- `3.visualize.py`: exports a VTF and multi-model PDB showing the ligand plus
  atomistic MAP pocket sidechains
- `4.compare_forcefield.py`: recomputes the imported ligand physical terms on the
  prepared snapshots and compares them to the live C++ node outputs

Run the full workflow:

```bash
source /Users/yinhan/Documents/upside2-md-test/.venv/bin/activate
source /Users/yinhan/Documents/upside2-md-test/source.sh
cd /Users/yinhan/Documents/upside2-md-test/example/17.ligand
python3 0.prepare.py
python3 1.run.py
python3 2.validate.py
python3 3.visualize.py
python3 4.compare_forcefield.py
```

or:

```bash
cd /Users/yinhan/Documents/upside2-md-test/example/17.ligand
./run.sh
```

Notes:
- this example minimizes shared-code changes by keeping ligand patching and
  visualization local to the example
- the smoke workflow applies example-local ligand-to-pocket distance restraints
  around the starting bound pose and runs with `--disable-recentering`, so this example is
  testing bound-pose refinement rather than unrestrained binding/unbinding
- pocket side chains are visualized as atomistic heavy-atom sidechains chosen from the
  ligand-conditioned MAP rotamer and placed with the live `affine_alignment` output
- `4.compare_forcefield.py` is an internal frozen-geometry consistency check; it does
  not replace comparison against OpenMM/CHARMM, but it verifies that the imported HDF5
  force-field data and the active C++ ligand nodes agree on the same snapshots
- the committed `charmm/` files are a minimal CHARMM-format subset for this example,
  not a full distributed CHARMM36m/CGenFF toppar package
- to point the example at a fuller local CHARMM/CGenFF install, set:
  - `UPSIDE_CHARMM_TOPPAR_DIR` to the local toppar root
  - optionally `UPSIDE_CHARMM_TOPOLOGY_FILES`, `UPSIDE_CHARMM_PARAMETER_FILES`,
    and `UPSIDE_CHARMM_STREAM_FILES` as `:`-separated relative or absolute paths
  - if those variables are unset, the workflow first uses `example/17.ligand/toppar`
    when present and otherwise falls back to the bundled example subset
- the example still expects the normal Upside Python dependencies used by the
  rest of the repo
