#!/bin/bash
# Minimal CG DOPC lipid test — builds tables, creates .up file, runs short minimization
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

source "${PROJECT_ROOT}/.venv/bin/activate"
source "${PROJECT_ROOT}/source.sh"

export UPSIDE_HOME="${PROJECT_ROOT}"
export PATH="${PROJECT_ROOT}/obj:$PATH"
export PYTHONPATH="${PROJECT_ROOT}/py:${PYTHONPATH:-}"

TEST_DIR="${SCRIPT_DIR}"
PDB_ID="test_cg_lipid"
RUNTIME_PDB="${TEST_DIR}/test_dopc_gly.pdb"
RUN_DIR="${TEST_DIR}/run"
MARTINI_H5="${RUN_DIR}/martini.h5"
UP_FILE="${RUN_DIR}/${PDB_ID}.stage_6.0.up"
CHECKPOINT_DIR="${RUN_DIR}/checkpoints"
UPSIDE_BIN="${PROJECT_ROOT}/obj/upside"

echo "=== CG Lipid Test Workflow ==="
echo "Project root: ${PROJECT_ROOT}"
echo "Test dir:     ${TEST_DIR}"
echo "Run dir:      ${RUN_DIR}"
echo ""

# --- Step 1: Build martini.h5 tables ---
echo "=== Step 1: Build martini tables ==="
mkdir -p "${RUN_DIR}"

export UPSIDE_MARTINI_FF_DIR="${PROJECT_ROOT}/parameters/dryMARTINI"

python3 -c "
import os, sys, numpy as np
from pathlib import Path
sys.path.insert(0, '${PROJECT_ROOT}/py')
from martini_build_tables import build_martini_tables
from martini_prepare_system_lib import _DOPC_ATOM_NAMES

# Parse DOPC reference bead positions from the test PDB
def parse_pdb_dopc(path):
    atoms = []
    with open(path, 'r') as f:
        for line in f:
            if not line.startswith(('ATOM','HETATM')):
                continue
            resname = line[17:21].strip().upper()
            if resname not in ('DOP', 'DOPC'):
                continue
            aname = line[12:16].strip().upper()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            atoms.append({'name': aname, 'x': x, 'y': y, 'z': z})
    return atoms

atoms = parse_pdb_dopc('${RUNTIME_PDB}')

# Group by residue (every 14 atoms)
n_per_lipid = 14
n_lipids = len(atoms) // n_per_lipid
if n_lipids == 0:
    raise RuntimeError('No DOPC lipids found in PDB')

# Use the first lipid's beads as reference (in Angstrom)
first_lipid = atoms[:n_per_lipid]
# Get COM of first lipid
com = np.mean([[a['x'], a['y'], a['z']] for a in first_lipid], axis=0)
# Reference bead positions relative to COM (in Angstrom, convert to nm for tables)
ref_bead_positions = np.array([[a['x']-com[0], a['y']-com[1], a['z']-com[2]] for a in first_lipid])
ref_bead_positions_nm = ref_bead_positions * 0.1  # Angstrom -> nm

print(f'Found {n_lipids} DOPC lipids, using first as reference')
print(f'Reference bead positions (nm): shape={ref_bead_positions_nm.shape}')

# DOPC bead types in ITP order
bead_types = ['Q0', 'Qa', 'Na', 'Na',
              'C1', 'C1', 'C3', 'C1', 'C1',
              'C1', 'C1', 'C3', 'C1', 'C1']

cg_lipid_config = {
    'ref_bead_positions_nm': ref_bead_positions_nm,
    'bead_types': bead_types,
}

ff_dir = Path(os.environ['UPSIDE_MARTINI_FF_DIR'])
dry_ff = ff_dir / 'dry_martini_v2.1.itp'
martinize = Path('${PROJECT_ROOT}/SC-training/') / 'martinize_dry.py'
sc_lib = Path('${PROJECT_ROOT}/SC-training/') / 'sidechain_orientation_library.h5'

print(f'Forcefield: {dry_ff}')
print(f'SC lib: {sc_lib}')

build_martini_tables(
    output_path='${MARTINI_H5}',
    dry_ff_path=dry_ff,
    martinize_path=martinize,
    sidechain_lib_path=sc_lib,
    forcefield_name='martini22',
    cg_lipid_config=cg_lipid_config,
)
print('Martini tables built successfully.')
"

echo ""

# --- Step 2: Convert stage (create .up file) ---
echo "=== Step 2: Convert stage ==="

export UPSIDE_RUNTIME_PDB_FILE="${RUNTIME_PDB}"
export UPSIDE_SIMULATION_STAGE="minimization"
export UPSIDE_MARTINI_FF_DIR="${PROJECT_ROOT}/parameters/dryMARTINI"

python3 -c "
import os, sys
sys.path.insert(0, '${PROJECT_ROOT}/py')
from martini_prepare_system_lib import convert_stage

convert_stage(
    pdb_id='${PDB_ID}',
    stage='minimization',
    run_dir='${RUN_DIR}',
)
print('Stage conversion complete.')
"

# --- Step 3: Inject particles table and CG lipid nodes ---
echo "=== Step 3: Inject tables ==="

python3 -c "
import sys
sys.path.insert(0, '${PROJECT_ROOT}/py')
from martini_prepare_system_lib import inject_particles_table, inject_cg_lipid_nodes
from pathlib import Path

up_file = Path('${UP_FILE}')
m5 = Path('${MARTINI_H5}')

inject_particles_table(up_file, m5)
print('Particles table injected.')

inject_cg_lipid_nodes(up_file, m5)
print('CG lipid nodes injected.')
"

# --- Step 4: Run short minimization ---
echo "=== Step 4: Run minimization (200 steps) ==="
mkdir -p "${CHECKPOINT_DIR}"

cat > "${RUN_DIR}/minimize.up" << 'UPEOF'
[global]
n_atom = 1028
temperature = 0.8647
time_step = 0.010
barostat = none
thermostat = langevin
thermostat_timescale = 5.0
seed = 42
max_iter = 200
energy_reporter_freq = 10
trajectory_reporter_freq = 0

[input]
file = ${PDB_ID}.stage_6.0.up
format = hdf5

[output]
checkpoint_dir = checkpoints
energy_file = energy.dat
trajectory_file = trajectory.xtc
UPEOF

echo "Run directory: ${RUN_DIR}"
echo "Done preparing test. To run minimization:"
echo "  cd ${RUN_DIR} && ${UPSIDE_BIN} minimize.up"
