# Hybrid Validation Notes

## Scope
This note records the Phase 7 validation checks for the hybrid dry-MARTINI + Upside integration implemented in this directory.

## What Was Validated

1. C++ build integrity
- Command:
```bash
cmake --build ../../obj -j4
```
- Result: build completed for `upside`, `upside_engine`, and `upside_calculation`.

2. Hybrid mapping schema and index consistency
- Generator command:
```bash
source ../../.venv/bin/activate
source ../../source.sh
python3 prepare_hybrid_system.py \
  --protein-pdb pdb/1rkl.pdb \
  --protein-cg-pdb pdb/1rkl.MARTINI.pdb \
  --bilayer-pdb pdb/bilayer.MARTINI.pdb \
  --output-dir outputs/hybrid_1rkl_test5
```
- Validation command:
```bash
source ../../.venv/bin/activate
python3 validate_hybrid_mapping.py \
  outputs/hybrid_1rkl_test5/hybrid_mapping.h5 \
  --n-atom 1263
```
- Result:
  - `n_atom=1263`
  - `n_protein=68`
  - `n_env=1195`
  - `n_bb=31`
  - `n_sc_rotamer_rows=37`

3. Runtime parser + coupling smoke test
- Created a small synthetic hybrid block on a copy of a runnable stage file:
  - `outputs/martini_test/checkpoints/bilayer.stage_7.0.hybrid_smoke.up`
- Ran a 1-step simulation:
```bash
source ../../source.sh
source ../../.venv/bin/activate
../../obj/upside outputs/martini_test/checkpoints/bilayer.stage_7.0.hybrid_smoke.up \
  --duration 1 \
  --frame-interval 1 \
  --temperature 0.8647 \
  --time-step 0.02 \
  --thermostat-timescale 4.0 \
  --thermostat-interval -1 \
  --seed 7090685331 \
  --integrator v \
  --disable-recentering
```
- Result:
  - Successful run completion.
  - Hybrid parser log confirmed:
    - `activation_stage=production`
    - `preprod_mode=rigid`
    - `n_bb=1`
    - `n_env=1019`

## Notes and Current Limits

1. Sidechain mapping is currently deterministic in prep output (`rotamer_probability=1.0` per proxy row). Full per-step probabilistic coupling from Upside rotamer populations is still future work.

2. Coupling-frame alignment (if enabled) currently uses a BB-frame rigid transform between consecutive frames (not full all-atom Kabsch).

3. Smoke test uses synthetic hybrid groups for parser/runtime coverage only; it is not a physical-quality validation trajectory.
