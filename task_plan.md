# Claude Code Plan: MARTINI Workflow Correction & NPT Logging

This plan outlines the steps to fix the multi-stage MARTINI simulation workflow, ensuring the NPT stage correctly utilizes the barostat and provides necessary debug logging for box stability.

---

## 1. Workflow Sequence Correction (run_sim_bilayer.sh)

The script `example/16.MARTINI/run_sim_bilayer.sh` must be updated to follow a strict sequential dependency where each stage uses the result of the previous stage as its starting structure.

**Required Sequence:**
1. **Prepare Input:** Generate the initial topology and coordinates.
2. **Minimization:** Relax the system to remove high-energy clashes.
3. **NPT Equilibration:** Allow the box dimensions to fluctuate to reach target pressure.
4. **NVT Production:** Run production at a fixed volume determined by the NPT stage.
5. **VTF Extraction:** Generate distinct trajectory files for both NPT and NVT stages.

### Current Workflow Analysis

**Stage 1: Prepare Input (Lines 83-100)**
- Runs `prepare_martini.py` with PDB ID
- Creates `.up` HDF5 file with system configuration
- Sets barostat configuration in `/input/barostat` group
- Output: `inputs/bilayer.up`

**Stage 2: Energy Minimization (Lines 105-141)**
- Duration: `MIN_STEPS` (default 500 steps)
- Uses softened potentials:
  - `UPSIDE_SOFTEN_LJ=1` with `UPSIDE_LJ_ALPHA=0.2`
  - `UPSIDE_SOFTEN_COULOMB=1` with `UPSIDE_SLATER_ALPHA=2.0`
- Runs with `--minimize` flag
- **Modifies input file in-place**
- Output: `checkpoints/bilayer.minimized.up` (copy of modified input)

**Stage 3: NPT Equilibration (Lines 146-189)**
- Duration: `NPT_STEPS` (default 2000 steps)
- NPT Configuration via Environment Variables:
  - `UPSIDE_NPT_ENABLE=1` (enables barostat)
  - `UPSIDE_NPT_TARGET_PXY=1.0` (lateral pressure target in bar)
  - `UPSIDE_NPT_TARGET_PZ=1.0` (normal pressure target in bar)
  - `UPSIDE_NPT_TAU=1.0` (barostat time constant in ps)
  - `UPSIDE_NPT_INTERVAL=10` (apply barostat every 10 steps)
- Uses `set_initial_position.py` to transfer last frame from minimization
- Output: `checkpoints/bilayer.npt.up`

**Stage 4: NVT Production (Lines 194-231)**
- Duration: `NVT_STEPS` (default 5000 steps)
- **NPT is disabled:** `UPSIDE_NPT_ENABLE=0`
- Uses fixed volume from NPT equilibration
- Uses `set_initial_position.py` to transfer last frame from NPT
- Output: `checkpoints/bilayer.nvt.up`

**Stage 5: VTF Generation (Lines 236-273)**
- Generates separate VTF files for each stage:
  - `bilayer.minimized.vtf`
  - `bilayer.npt.vtf`
  - `bilayer.nvt.vtf`
- Uses `extract_martini_vtf.py` for trajectory extraction

---

## 2. Implementation Tasks

### Task A: Fix Stage Logic & Continuity

#### A1. NPT Activation Verification ✅ VERIFIED
**Finding:** NPT is configured via environment variables, NOT command-line flags.

**How NPT Works:**
1. Environment variables (e.g., `UPSIDE_NPT_ENABLE=1`) are read by `prepare_martini.py`
2. Settings are stored in HDF5 file at `/input/barostat` group
3. C++ code (`box.cpp:register_barostat_for_engine()`) reads settings at runtime
4. No command-line flags like `--npt` are needed

**Current Implementation:**
- `prepare_martini.py` (lines 1104-1122): Creates `/input/barostat` group with attributes
- `box.cpp` (lines 147-189): Reads barostat settings from HDF5
- `main.cpp` (lines 1097-1103): Calls `maybe_apply_barostat()` during integration

**Action Required:** None - NPT activation mechanism is correct.

#### A2. Sequential Continuity ✅ VERIFIED
**Finding:** Stage transitions are correctly implemented.

**Current Mechanism:**
```bash
# Stage 2 → 3
cp -f "$MINIMIZED_FILE" "$NPT_FILE"
python3 set_initial_position.py "$MINIMIZED_FILE" "$NPT_FILE"

# Stage 3 → 4
cp -f "$NPT_FILE" "$NVT_FILE"
python3 set_initial_position.py "$NPT_FILE" "$NVT_FILE"
```

**What `set_initial_position.py` Does:**
- Reads last frame from `/output/pos[-1]` of source file
- Writes to `/input/pos` of destination file
- Preserves all other configuration (topology, potentials, barostat settings)

**Action Required:** None - continuity is correct.

#### A3. VTF Separation ✅ VERIFIED
**Finding:** VTF files are already generated separately for each stage.

**Current Implementation:**
- Stage 5 calls `extract_martini_vtf.py` three times with different input files
- Each stage produces a unique VTF file
- Files are named: `bilayer.minimized.vtf`, `bilayer.npt.vtf`, `bilayer.nvt.vtf`

**Action Required:** None - VTF separation is correct.

### Task B: Source Code Implementation (NPT Debugging)

#### B1. Current Logging Status

**Console Output (Already Implemented):**
The barostat already prints debug information to stdout when `debug=true`:

```cpp
// Registration message (box.cpp:186)
[NPT] Barostat registered: box 50.00 x 50.00 x 50.00, target Pxy=1.000e+00 Pz=1.000e+00, tau=1.00, interval=10

// Per-application message (box.cpp:299)
 [NPT] t 0.100 scale_xy 1.0023 scale_z 0.9987 | Pxy 1.05e+00 tgt 1.00e+00, Pz 9.87e-01 tgt 1.00e+00 | box 50.12 50.12 49.94
```

**HDF5 Output (NOT Implemented):**
Currently, box dimensions and pressure are **NOT logged to HDF5 output**, only printed to console.

#### B2. Missing HDF5 Datasets

The following NPT-related data is calculated but not logged to HDF5:

1. **Box Dimensions Over Time**
   - Current values: `box_x`, `box_y`, `box_z` (stored in `BarostatState`)
   - Should be logged as `/output/box` dataset with shape `[n_frames, 3]`

2. **Instantaneous Pressure**
   - Lateral pressure: `pxy_inst` (calculated at box.cpp:222)
   - Normal pressure: `pz_inst` (calculated at box.cpp:223)
   - Should be logged as `/output/pressure` dataset with shape `[n_frames, 2]`

3. **Pressure Tensor Components**
   - Components: `pxx`, `pyy`, `pzz` (calculated at box.cpp:220)
   - Should be logged as `/output/pressure_tensor` dataset with shape `[n_frames, 3]`

4. **Scale Factors**
   - Lateral scale: `scale_xy` (calculated at box.cpp:228-257)
   - Normal scale: `scale_z` (calculated at box.cpp:228-257)
   - Should be logged as `/output/scale_factors` dataset with shape `[n_frames, 2]`

5. **Volume**
   - Calculated at box.cpp:212 but not stored
   - Should be logged as `/output/volume` dataset with shape `[n_frames]`

#### B3. Implementation Approach

**Option A: Add Logging in main.cpp (Recommended)**

**Location:** After temperature logger (main.cpp:927-933)

**Steps:**
1. Add accessor functions in `box.h`/`box.cpp` to retrieve:
   - Current box dimensions
   - Last calculated pressure values
   - Last scale factors
   - Current volume

2. Store pressure and scale values in `BarostatState` (box.h:54-64):
   ```cpp
   struct BarostatState {
       // ... existing fields ...
       float last_pxy_inst;
       float last_pz_inst;
       float last_scale_xy;
       float last_scale_z;
   };
   ```

3. Update `maybe_apply_barostat()` to store these values in `BarostatState`

4. Add loggers in main.cpp:
   ```cpp
   // After line 933
   if (simulation_box::npt::is_enabled(sys.engine)) {
       sys.logger.add_logger<float>("box", {3},
           [&sys]() {
               float bx, by, bz;
               simulation_box::npt::get_current_box(sys.engine, bx, by, bz);
               return vector<float>{bx, by, bz};
           });

       sys.logger.add_logger<float>("pressure", {2},
           [&sys]() {
               float pxy, pz;
               simulation_box::npt::get_pressure(sys.engine, pxy, pz);
               return vector<float>{pxy, pz};
           });

       sys.logger.add_logger<float>("volume", {},
           [&sys]() {
               return simulation_box::npt::get_volume(sys.engine);
           });
   }
   ```

**Option B: Add Logging in box.cpp**

**Location:** Inside `maybe_apply_barostat()` after pressure calculation

**Pros:**
- Direct access to all calculated values
- No need for accessor functions

**Cons:**
- Requires passing logger to `maybe_apply_barostat()`
- More invasive changes to function signature
- Harder to maintain separation of concerns

**Recommendation:** Use Option A for cleaner architecture.

#### B4. Detailed Implementation Plan

**File 1: `src/box.h`**
- Add fields to `BarostatState` struct (lines 54-64):
  - `float last_pxy_inst`
  - `float last_pz_inst`
  - `float last_scale_xy`
  - `float last_scale_z`
- Add accessor function declarations in `simulation_box::npt` namespace:
  - `void get_pressure(const DerivEngine& engine, float& pxy, float& pz)`
  - `float get_volume(const DerivEngine& engine)`

**File 2: `src/box.cpp`**
- Update `maybe_apply_barostat()` (lines 191-305):
  - Store `pxy_inst` in `st.last_pxy_inst` (after line 222)
  - Store `pz_inst` in `st.last_pz_inst` (after line 223)
  - Store `scale_xy` in `st.last_scale_xy` (after line 257)
  - Store `scale_z` in `st.last_scale_z` (after line 257)
- Implement accessor functions (after line 325):
  ```cpp
  void get_pressure(const DerivEngine& engine, float& pxy, float& pz) {
      auto& st = engine.get_computation<BarostatState&>("barostat_state");
      pxy = st.last_pxy_inst;
      pz = st.last_pz_inst;
  }

  float get_volume(const DerivEngine& engine) {
      float bx, by, bz;
      get_current_box(engine, bx, by, bz);
      return bx * by * bz;
  }
  ```

**File 3: `src/main.cpp`**
- Add NPT loggers after temperature logger (after line 933):
  ```cpp
  // Add NPT-specific loggers if barostat is enabled
  if (simulation_box::npt::is_enabled(sys.engine)) {
      sys.logger.add_logger<float>("box", {3},
          [&sys]() {
              float bx, by, bz;
              simulation_box::npt::get_current_box(sys.engine, bx, by, bz);
              return vector<float>{bx, by, bz};
          });

      sys.logger.add_logger<float>("pressure", {2},
          [&sys]() {
              float pxy, pz;
              simulation_box::npt::get_pressure(sys.engine, pxy, pz);
              return vector<float>{pxy, pz};
          });

      sys.logger.add_logger<float>("volume", {},
          [&sys]() {
              return simulation_box::npt::get_volume(sys.engine);
          });
  }
  ```

---

## 3. Verification & Quality Control

### V1. Build Verification
```bash
cd /Users/yinhan/Documents/upside2-md-dev
./install_M1.sh
```
**Expected:** Build completes without errors

### V2. NPT Activation Check
```bash
cd example/16.MARTINI
./run_sim_bilayer.sh
```
**Expected:** Console output shows `[NPT]` messages with changing box dimensions

### V3. HDF5 Output Verification
```bash
# Check that new datasets exist
h5ls -r checkpoints/bilayer.npt.up | grep -E "box|pressure|volume"
```
**Expected Output:**
```
/output/box                      Dataset {n_frames, 3}
/output/pressure                 Dataset {n_frames, 2}
/output/volume                   Dataset {n_frames}
```

### V4. Box Dimension Stability Check
```python
import h5py
import numpy as np

with h5py.File('checkpoints/bilayer.npt.up', 'r') as f:
    box = f['/output/box'][:]
    pressure = f['/output/pressure'][:]

    # Check that box dimensions change during NPT
    box_std = np.std(box, axis=0)
    print(f"Box dimension std dev: {box_std}")
    assert np.all(box_std > 0.01), "Box dimensions should fluctuate during NPT"

    # Check that pressure converges to target
    target_pxy = 1.0
    target_pz = 1.0
    final_pxy = np.mean(pressure[-100:, 0])
    final_pz = np.mean(pressure[-100:, 1])
    print(f"Final pressure: Pxy={final_pxy:.3f}, Pz={final_pz:.3f}")
    assert abs(final_pxy - target_pxy) < 0.1, "Lateral pressure should converge"
    assert abs(final_pz - target_pz) < 0.1, "Normal pressure should converge"
```

### V5. VTF Extraction Verification
```bash
# Check that VTF files use correct box dimensions
grep "pbc" outputs/bilayer.npt.vtf | head -1
grep "pbc" outputs/bilayer.nvt.vtf | head -1
```
**Expected:** NPT and NVT VTF files should have different box dimensions if equilibration occurred

### V6. Dry Run Test
```bash
# Run short simulation to test logging
cd example/16.MARTINI
MIN_STEPS=10 NPT_STEPS=100 NVT_STEPS=10 ./run_sim_bilayer.sh
```
**Expected:**
- Simulation completes in <1 minute
- Console shows `[NPT]` messages
- HDF5 files contain new datasets
- No formatting errors in logs

---

## 4. Critical Files to Modify

### Source Code Files
1. **src/box.h** (lines 54-64, after line 325)
   - Add fields to `BarostatState` struct
   - Add accessor function declarations

2. **src/box.cpp** (lines 191-305, after line 325)
   - Store pressure and scale values in `BarostatState`
   - Implement accessor functions

3. **src/main.cpp** (after line 933)
   - Add NPT loggers for box, pressure, and volume

### Workflow Files (No Changes Needed)
4. **example/16.MARTINI/run_sim_bilayer.sh** ✅ Already correct
5. **example/16.MARTINI/set_initial_position.py** ✅ Already correct
6. **example/16.MARTINI/extract_martini_vtf.py** ✅ Already correct

---

## 5. Known Issues & Limitations

### Issue 1: VTF Extraction Fallback Logic
**Problem:** `extract_martini_vtf.py` has complex fallback logic for finding box dimensions:
1. `/output/box` dataset (will now exist after our changes)
2. Log file parsing (brittle, depends on format)
3. `/equilibrated_box` group attributes
4. `/input/potential/martini_potential` attributes
5. PDB file CRYST1 record
6. Default fallback: 50.0 x 50.0 x 50.0 Å

**Impact:** After implementing HDF5 logging, VTF extraction will use `/output/box` dataset (priority 1), making it more robust.

### Issue 2: NPT Configuration Persistence
**Problem:** NPT settings are baked into HDF5 file during preparation. Changing NPT settings requires re-running `prepare_martini.py`.

**Workaround:** Use environment variables before running `prepare_martini.py`:
```bash
export UPSIDE_NPT_TARGET_PXY=2.0
export UPSIDE_NPT_TARGET_PZ=1.5
python3 prepare_martini.py ...
```

### Issue 3: Equilibrium Detection
**Problem:** No automatic detection of when NPT equilibration has converged.

**Workaround:** Manually inspect `/output/pressure` dataset to verify convergence before starting NVT production.

---

## 6. Execution Summary

### What Works (No Changes Needed)
✅ NPT activation via environment variables
✅ Sequential stage transitions with `set_initial_position.py`
✅ Separate VTF file generation for each stage
✅ Console debug output for NPT

### What Needs Implementation
❌ HDF5 logging of box dimensions
❌ HDF5 logging of pressure values
❌ HDF5 logging of volume
❌ Accessor functions for NPT state

### Implementation Effort
- **Files to modify:** 3 (box.h, box.cpp, main.cpp)
- **Lines to add:** ~50 lines total
- **Complexity:** Low (adding loggers to existing infrastructure)
- **Risk:** Low (no changes to core simulation logic)

---

## 7. Post-Implementation Benefits

1. **Robust VTF Extraction:** Box dimensions will be read from HDF5 instead of log file parsing
2. **NPT Analysis:** Pressure and volume trajectories available for convergence analysis
3. **Debugging:** Easy to verify that NPT is working correctly
4. **Reproducibility:** All NPT data stored in HDF5 for future reference
5. **Visualization:** Can plot box dimensions and pressure over time
