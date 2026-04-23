# Fix: Preserve Multi-Chain Topology For 1AFO Prep

## Project Goal
- Make the `1AFO` preparation path preserve its two peptide chains instead of collapsing them into one continuous chain in the generated Upside stage files.
- Ensure both metadata and generated stage-7 backbone terms reflect the actual chain break between chain `A` and chain `B`.

## Architecture & Key Decisions
- Keep the fix in the shared hybrid metadata and stage-7 injection path rather than in the thin `1afo` wrappers.
- Write explicit chain-break metadata from the AA-backbone residue order into the hybrid metadata HDF5.
- Copy that chain-break metadata into prepared stage files along with the existing hybrid groups.
- Make backbone-node injection read the chain-break metadata and configure Upside with the correct `n_chains` / `chain_starts` instead of forcing a single chain.
- Also persist per-residue chain IDs in the backbone map so downstream exports can recover the correct protein chain labels.

## Execution Phases
- [x] Phase 1: Update the trackers for the multi-chain prep bug and inspect the current metadata / backbone-node injection path.
- [x] Phase 2: Patch the metadata writer, stage metadata injection, and backbone-node injection so multi-chain proteins preserve chain breaks.
- [x] Phase 3: Verify the `1AFO` metadata and a reduced generated stage file show the expected two-chain topology and no peptide bond across the chain boundary.

## Known Errors / Blockers
- No blocker.

## Review
- The shared hybrid metadata writer now preserves per-residue chain IDs and writes `/input/chain_break` for multi-chain proteins.
- The stage metadata handoff now copies `chain_break` into prepared stage files.
- Backbone-node injection no longer forces `n_chains = 1`; it reads the stored chain break and configures Upside with the correct chain starts.
- `martini_extract_vtf.py` now preserves AA-backbone chain IDs and skips cross-chain backbone bonds in the exported VTF topology.
- Verified with:
  - `PYTHONPYCACHEPREFIX=/tmp python3 -m py_compile ../../py/martini_prepare_system_lib.py ../../py/martini_extract_vtf.py`
  - `bash -n run_sim_1rkl.sh`
  - metadata round-trip on `/tmp/1afo_chainbreak_metadata.h5`
  - reduced local workflow run:
    - `RUN_DIR=outputs/1afo_chainbreak_check ... bash run_sim_1afo.sh`
  - generated production file checks:
    - `/input/chain_break/chain_first_residue = [36]`
    - `/input/chain_break/chain_counts = [1, 1]`
    - no `Distance3D` bond between the last chain-`A` C atom and the first chain-`B` N atom
  - generated VTF checks:
    - protein entries appear under both chain `A` and chain `B`.
