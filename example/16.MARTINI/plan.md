# Workflow Fix: 1AFO Sidechain Injection Residue Collapse

## Project Goal
- Fix the `1afo_outlipid` stage-preparation failure:
  - `ValueError: Missing or inconsistent /input/sequence for AA-backbone sidechain injection: expected 36 residues`
- Preserve the existing `1rkl` behavior while making multi-chain AA-backbone workflows like `1AFO` generate consistent hybrid metadata.

## Architecture & Key Decisions
- Fix the problem in the shared metadata-generation path, not in the `1afo` wrapper scripts.
- The root cause is a residue-identity collapse in `hybrid_bb_map/bb_residue_index`:
  - sequence extraction is chain-aware,
  - backbone-map residue IDs currently use raw `resseq` only,
  - multi-chain proteins with repeated residue numbers therefore collapse distinct residues.
- The fix should make `bb_residue_index` unique per backbone residue in residue-order space.
- Keep raw chain/resseq information only as descriptive metadata/comments; do not rely on raw `resseq` as the unique hybrid residue key.

## Execution Phases
- [x] Phase 1: Confirm the multi-chain residue-index collapse and identify the shared metadata write path to change.
- [x] Phase 2: Patch the hybrid metadata generation so backbone residue IDs remain unique across chains.
- [x] Phase 3: Re-run the failing `1AFO` prep/injection path and verify the sidechain injection succeeds.

## Known Errors / Blockers
- No blocker.

## Review
- Root cause confirmed in the shared hybrid metadata path:
  - `extract_backbone_sequence(...)` is chain-aware and returns `72` residues for `1AFO`,
  - `collect_aa_backbone_map(...)` previously stored raw PDB `resseq` into `bb_residue_index`,
  - because chains `A` and `B` both use residue numbers `66..101`, the written metadata collapsed to only `36` unique backbone residue IDs,
  - `inject_stage7_sc_table_nodes(...)` therefore saw `72` sequence entries but only `36` affine residues and aborted.
- Implemented fix:
  - `py/martini_prepare_system_lib.py::collect_aa_backbone_map(...)` now writes a residue-order `bb_residue_index` that is unique per backbone residue across chains,
  - `write_backbone_metadata_h5(...)` now writes that unique index into `hybrid_bb_map/bb_residue_index`, with backward-compatible fallback to the old `bb_resseq` field if needed.
- Verification:
  - helper check on `example/16.MARTINI/pdb/1AFO.pdb`:
    - `len(entries) = 72`
    - `len(sequence) = 72`
    - `unique bb_residue_index = 72`
  - metadata write/read check on `/tmp/1afo_test_backbone_metadata.h5`:
    - `sequence len = 72`
    - `unique bb_residue_index = 72`
  - reduced local workflow smoke:
    - `RUN_DIR=example/16.MARTINI/outputs/1afo_inject_fix_check`
    - `AUTO_CONTINUE_FROM_PREVIOUS_RUN=0`
    - `DISABLE_1AFO_AABB_AUTO_CONTINUE=1`
    - minimal stage lengths set to `1`
    - the run advanced past the previous failure site:
      - stage-0 packing succeeded,
      - stage conversion wrote `test.input.up`,
      - `inject-stage7-sc` no longer raised `expected 36 residues`.
