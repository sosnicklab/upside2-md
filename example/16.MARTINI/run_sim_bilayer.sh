#!/bin/bash
source ../../source.sh
source ../../.venv/bin/activate

# MARTINI 3.0 Bilayer Simulation Workflow (CHARMM-GUI Protocol)
# Stages: Prepare -> Minimization -> NPT Equilibration (Berendsen) -> NPT Production (Parrinello-Rahman) -> VTF Generation
# Each stage uses the output checkpoint from the previous stage as its starting point.
#
# ENVIRONMENT VARIABLES:
#   PDB_ID                - PDB identifier (default: bilayer)
#   UPSIDE_NPT_TARGET_PXY - Target lateral pressure (default: 0.000020659 = 1 bar)
#   UPSIDE_NPT_TARGET_PZ  - Target normal pressure (default: 0.000020659 = 1 bar)
#   UPSIDE_NPT_TAU        - Barostat time constant (default: 1.0)
#   UPSIDE_NPT_INTERVAL   - Steps between barostat applications (default: 10)
#   UPSIDE_BAROSTAT_TYPE  - Barostat type: 0=Berendsen, 1=Parrinello-Rahman

set -e  # Exit on any error

# =============================================================================
# USER CONFIGURATION
# =============================================================================
PDB_ID="${PDB_ID:-bilayer}"

# Directories
INPUTS_DIR="inputs"
OUTPUTS_DIR="outputs"
RUN_DIR="outputs/martini_test"
CHECKPOINT_DIR="${RUN_DIR}/checkpoints"

# Filenames
INPUT_FILE="${INPUTS_DIR}/${PDB_ID}.up"
MINIMIZED_FILE="${CHECKPOINT_DIR}/${PDB_ID}.minimized.up"
NPT_EQUIL_FILE="${CHECKPOINT_DIR}/${PDB_ID}.npt_equil.up"
NPT_PROD_FILE="${CHECKPOINT_DIR}/${PDB_ID}.npt_prod.up"
VTF_FILE="${RUN_DIR}/${PDB_ID}.vtf"
LOG_DIR="${RUN_DIR}/logs"

# Simulation parameters
TEMPERATURE=0.8
TIME_STEP=0.1
THERMOSTAT_TIMESCALE=0.135
THERMOSTAT_INTERVAL="${THERMOSTAT_INTERVAL:--1}"
SEED=12345

# Stage durations (in MD steps)
MIN_STEPS="${MIN_STEPS:-500}"
NPT_EQUIL_STEPS="${NPT_EQUIL_STEPS:-2000}"
NPT_PROD_STEPS="${NPT_PROD_STEPS:-5000}"
FRAME_INTERVAL="${FRAME_INTERVAL:-20}"

# Softening parameters for minimization
export UPSIDE_SOFTEN_LJ=${UPSIDE_SOFTEN_LJ:-1}
export UPSIDE_LJ_ALPHA=${UPSIDE_LJ_ALPHA:-0.2}
export UPSIDE_SOFTEN_COULOMB=${UPSIDE_SOFTEN_COULOMB:-1}
export UPSIDE_SLATER_ALPHA=${UPSIDE_SLATER_ALPHA:-2.0}

# =============================================================================
# VALIDATION
# =============================================================================
if [ -z "$UPSIDE_HOME" ]; then
    echo "ERROR: UPSIDE_HOME environment variable is not set!"
    exit 1
fi

UPSIDE_EXECUTABLE="${UPSIDE_HOME}/obj/upside"
if [ ! -f "$UPSIDE_EXECUTABLE" ]; then
    echo "ERROR: UPSIDE executable not found: $UPSIDE_EXECUTABLE"
    exit 1
fi

# Create directories
mkdir -p "$INPUTS_DIR" "$OUTPUTS_DIR" "$RUN_DIR" "$CHECKPOINT_DIR" "$LOG_DIR"

echo "=== MARTINI 3.0 Bilayer Workflow (CHARMM-GUI Protocol) ==="
echo "PDB ID: $PDB_ID"
echo "Stages: Prepare -> Minimization -> NPT Equilibration (Berendsen) -> NPT Production (Parrinello-Rahman) -> VTF"
echo "  NPT Equilibration:  $NPT_EQUIL_STEPS steps"
echo "  NPT Production:     $NPT_PROD_STEPS steps"
echo

# =============================================================================
# STAGE 1: PREPARE INPUT FILES
# =============================================================================
echo "=== Stage 1: Preparing Input Files ==="
source ../../source.sh
source ../../.venv/bin/activate

# Prepare with NPT enabled (so barostat configuration is stored in input file)
export UPSIDE_OVERWRITE_SPLINES=${UPSIDE_OVERWRITE_SPLINES:-1}
export UPSIDE_NPT_ENABLE=${UPSIDE_NPT_ENABLE:-1}
python3 prepare_martini.py "$PDB_ID"

PREPARED_FILE="${RUN_DIR}/test.input.up"
if [ ! -f "$PREPARED_FILE" ]; then
    echo "ERROR: Input preparation failed - file not found: $PREPARED_FILE"
    exit 1
fi

# Copy to normalized location
cp -f "$PREPARED_FILE" "$INPUT_FILE"
echo "Input prepared: $INPUT_FILE ($(du -h "$INPUT_FILE" | cut -f1))"
echo

# =============================================================================
# STAGE 2: ENERGY MINIMIZATION
# =============================================================================
echo "=== Stage 2: Energy Minimization ==="
echo "Running gradient descent minimization with softened potentials"
echo "Input:  $INPUT_FILE"
echo "Output: $MINIMIZED_FILE"

CMD_MIN=(
    "$UPSIDE_EXECUTABLE"
    "$INPUT_FILE"
    # Required arguments even for minimization (executable parses all args before checking --minimize)
    "--duration" "0"          # 0 duration means no MD steps after minimization
    "--frame-interval" "1"    # Required by executable
    "--temperature" "$TEMPERATURE"  # Required by executable
    "--time-step" "$TIME_STEP"      # Required by executable
    "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"  # Required by executable
    "--thermostat-interval" "$THERMOSTAT_INTERVAL"    # Required by executable
    "--seed" "$SEED"          # Required by executable
    "--integrator" "vel_verlet"    # Required by executable
    "--disable-recentering"
    "--minimize"
    "--min-max-iter" "1000"
    "--min-energy-tol" "1e-6"
    "--min-force-tol" "1e-3"
    "--min-step" "0.01"
)

START_TIME=$(date +%s)
if "${CMD_MIN[@]}" 2>&1 | tee "${LOG_DIR}/minimization.log"; then
    END_TIME=$(date +%s)
    echo "Minimization completed in $((END_TIME - START_TIME)) seconds"
else
    echo "ERROR: Minimization failed!"
    exit 1
fi

# Save minimized checkpoint (input file was modified in-place)
cp -f "$INPUT_FILE" "$MINIMIZED_FILE"
echo "Minimized checkpoint: $MINIMIZED_FILE"
echo

# =============================================================================
# STAGE 3: NPT EQUILIBRATION (Berendsen barostat)
# =============================================================================
echo "=== Stage 3: NPT Equilibration (Berendsen) ==="
echo "Running $NPT_EQUIL_STEPS steps with Berendsen barostat"
echo "Input:  $MINIMIZED_FILE"
echo "Output: $NPT_EQUIL_FILE"

# Disable softened potentials for NPT equilibration (use hard particles)
# NOTE: We need to modify the HDF5 file attributes directly, not just environment variables
# because Upside reads these parameters from the HDF5 file, not from environment variables
export UPSIDE_SOFTEN_LJ=0
export UPSIDE_SOFTEN_COULOMB=0

# Enable NPT with Berendsen barostat for equilibration
export UPSIDE_NPT_ENABLE=1
export UPSIDE_BAROSTAT_TYPE=0  # Berendsen
# Pressure in Upside units: 1 bar = 0.000020659 E_up/Angstrom^3
export UPSIDE_NPT_TARGET_PXY=${UPSIDE_NPT_TARGET_PXY:-0.000020659}
export UPSIDE_NPT_TARGET_PZ=${UPSIDE_NPT_TARGET_PZ:-0.000020659}
export UPSIDE_NPT_TAU=${UPSIDE_NPT_TAU:-1.0}
export UPSIDE_NPT_INTERVAL=${UPSIDE_NPT_INTERVAL:-10}

echo "NPT settings: Pxy=${UPSIDE_NPT_TARGET_PXY} Pz=${UPSIDE_NPT_TARGET_PZ} tau=${UPSIDE_NPT_TAU} interval=${UPSIDE_NPT_INTERVAL}"
echo "Barostat type: Berendsen"
echo "Using hard particles (disabling softened potentials)"

# Work from minimized checkpoint
cp -f "$MINIMIZED_FILE" "$NPT_EQUIL_FILE"

# Disable softened potentials in the HDF5 file
python3 - <<END
import h5py

with h5py.File("$NPT_EQUIL_FILE", 'r+') as f:
    if '/input/potential/martini_potential' in f:
        grp = f['/input/potential/martini_potential']
        # Disable Coulomb softening (Slater softening)
        grp.attrs['coulomb_soften'] = 0
        # Disable LJ softening (soft-core LJ)
        grp.attrs['lj_soften'] = 0
        print("Disabled softened potentials in HDF5 file:")
        print(f"  coulomb_soften = {grp.attrs['coulomb_soften']}")
        print(f"  lj_soften = {grp.attrs['lj_soften']}")
END

# Set initial position in NPT file to last frame from minimization
python3 set_initial_position.py "$MINIMIZED_FILE" "$NPT_EQUIL_FILE"

CMD_NPT_EQUIL=(
    "$UPSIDE_EXECUTABLE"
    "$NPT_EQUIL_FILE"
    "--duration" "$NPT_EQUIL_STEPS"
    "--frame-interval" "$FRAME_INTERVAL"
    "--temperature" "$TEMPERATURE"
    "--time-step" "$TIME_STEP"
    "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
    "--thermostat-interval" "$THERMOSTAT_INTERVAL"
    "--seed" "$SEED"
    "--integrator" "vel_verlet"
    "--disable-recentering"
)

START_TIME=$(date +%s)
if "${CMD_NPT_EQUIL[@]}" 2>&1 | tee "${LOG_DIR}/npt_equilibration.log"; then
    END_TIME=$(date +%s)
    echo "NPT equilibration completed in $((END_TIME - START_TIME)) seconds"
else
    echo "ERROR: NPT equilibration failed!"
    exit 1
fi

echo "NPT equilibration checkpoint: $NPT_EQUIL_FILE"
echo

# =============================================================================
# STAGE 4: NPT PRODUCTION (Parrinello-Rahman barostat)
# =============================================================================
echo "=== Stage 4: NPT Production (Parrinello-Rahman) ==="
echo "Running $NPT_PROD_STEPS steps with Parrinello-Rahman barostat"
echo "Input:  $NPT_EQUIL_FILE"
echo "Output: $NPT_PROD_FILE"

# Switch to Parrinello-Rahman barostat for production
export UPSIDE_BAROSTAT_TYPE=1  # Parrinello-Rahman

echo "NPT settings: Pxy=${UPSIDE_NPT_TARGET_PXY} Pz=${UPSIDE_NPT_TARGET_PZ} tau=${UPSIDE_NPT_TAU} interval=${UPSIDE_NPT_INTERVAL}"
echo "Barostat type: Parrinello-Rahman"
echo "Using hard particles (softened potentials already disabled)"

# Work from equilibrated checkpoint
cp -f "$NPT_EQUIL_FILE" "$NPT_PROD_FILE"
# Update barostat type to Parrinello-Rahman AND update box dimensions from equilibrated state
python3 - <<END
import h5py
import numpy as np

with h5py.File("$NPT_PROD_FILE", 'r+') as f:
    # Update barostat type to Parrinello-Rahman
    if '/input/barostat' in f:
        f['/input/barostat'].attrs['type'] = 1  # Parrinello-Rahman

    # Ensure softened potentials are still disabled in production file
    if '/input/potential/martini_potential' in f:
        grp = f['/input/potential/martini_potential']
        # Double-check Coulomb softening is disabled
        grp.attrs['coulomb_soften'] = 0
        # Double-check LJ softening is disabled
        grp.attrs['lj_soften'] = 0
        print("Verified softened potentials are disabled in production file:")
        print(f"  coulomb_soften = {grp.attrs['coulomb_soften']}")
        print(f"  lj_soften = {grp.attrs['lj_soften']}")

    # Get equilibrated box dimensions from last frame of output
    if '/output/box' in f:
        last_box = f['/output/box'][-1]  # [x, y, z]
        print(f"Updating box dimensions from equilibrated state: {last_box}")

        # Update box dimensions in martini_potential attributes
        if '/input/potential/martini_potential' in f:
            grp = f['/input/potential/martini_potential']
            grp.attrs['x_len'] = float(last_box[0])
            grp.attrs['y_len'] = float(last_box[1])
            grp.attrs['z_len'] = float(last_box[2])
            print(f"Updated martini_potential: x_len={last_box[0]:.6f}, y_len={last_box[1]:.6f}, z_len={last_box[2]:.6f}")
END
# Set initial position in production file to last frame from equilibration
python3 set_initial_position.py "$NPT_EQUIL_FILE" "$NPT_PROD_FILE"

CMD_NPT_PROD=(
    "$UPSIDE_EXECUTABLE"
    "$NPT_PROD_FILE"
    "--duration" "$NPT_PROD_STEPS"
    "--frame-interval" "$FRAME_INTERVAL"
    "--temperature" "$TEMPERATURE"
    "--time-step" "$TIME_STEP"
    "--thermostat-timescale" "$THERMOSTAT_TIMESCALE"
    "--thermostat-interval" "$THERMOSTAT_INTERVAL"
    "--seed" "$SEED"
    "--integrator" "vel_verlet"
    "--disable-recentering"
)

START_TIME=$(date +%s)
if "${CMD_NPT_PROD[@]}" 2>&1 | tee "${LOG_DIR}/npt_production.log"; then
    END_TIME=$(date +%s)
    echo "NPT production completed in $((END_TIME - START_TIME)) seconds"
else
    echo "ERROR: NPT production failed!"
    exit 1
fi

echo "NPT production checkpoint: $NPT_PROD_FILE"
echo

# =============================================================================
# STAGE 5: VTF GENERATION (ALL STAGES)
# =============================================================================
echo "=== Stage 5: Generating VTF Files for All Stages ==="

# Minimization VTF
MIN_VTF_FILE="${RUN_DIR}/${PDB_ID}.minimized.vtf"
echo "Extracting minimization trajectory from: $MINIMIZED_FILE"
echo "Output VTF: $MIN_VTF_FILE"
if python3 extract_martini_vtf.py "$MINIMIZED_FILE" "$MIN_VTF_FILE" "$MINIMIZED_FILE" "$PDB_ID"; then
    MIN_VTF_SIZE=$(du -h "$MIN_VTF_FILE" | cut -f1)
    echo "Minimization VTF generated: $MIN_VTF_FILE ($MIN_VTF_SIZE)"
else
    echo "ERROR: Minimization VTF generation failed!"
    exit 1
fi

# NPT Equilibration VTF
NPT_EQUIL_VTF_FILE="${RUN_DIR}/${PDB_ID}.npt_equil.vtf"
echo "Extracting NPT equilibration trajectory from: $NPT_EQUIL_FILE"
echo "Output VTF: $NPT_EQUIL_VTF_FILE"
if python3 extract_martini_vtf.py "$NPT_EQUIL_FILE" "$NPT_EQUIL_VTF_FILE" "$NPT_EQUIL_FILE" "$PDB_ID"; then
    NPT_EQUIL_VTF_SIZE=$(du -h "$NPT_EQUIL_VTF_FILE" | cut -f1)
    echo "NPT equilibration VTF generated: $NPT_EQUIL_VTF_FILE ($NPT_EQUIL_VTF_SIZE)"
else
    echo "ERROR: NPT equilibration VTF generation failed!"
    exit 1
fi

# NPT Production VTF
NPT_PROD_VTF_FILE="${RUN_DIR}/${PDB_ID}.npt_prod.vtf"
echo "Extracting NPT production trajectory from: $NPT_PROD_FILE"
echo "Output VTF: $NPT_PROD_VTF_FILE"
if python3 extract_martini_vtf.py "$NPT_PROD_FILE" "$NPT_PROD_VTF_FILE" "$NPT_PROD_FILE" "$PDB_ID"; then
    NPT_PROD_VTF_SIZE=$(du -h "$NPT_PROD_VTF_FILE" | cut -f1)
    echo "NPT production VTF generated: $NPT_PROD_VTF_FILE ($NPT_PROD_VTF_SIZE)"
else
    echo "ERROR: NPT production VTF generation failed!"
    exit 1
fi

echo
echo "=== Workflow Complete ==="
echo "Checkpoints:"
echo "  Prepared:          $INPUT_FILE"
echo "  Minimized:         $MINIMIZED_FILE"
echo "  NPT Equilibration: $NPT_EQUIL_FILE"
echo "  NPT Production:    $NPT_PROD_FILE"
echo "Trajectories:"
echo "  Minimization:      $MIN_VTF_FILE"
echo "  NPT Equilibration: $NPT_EQUIL_VTF_FILE"
echo "  NPT Production:    $NPT_PROD_VTF_FILE"
echo
echo "To visualize: vmd <vtf_file>"
echo "Done."
