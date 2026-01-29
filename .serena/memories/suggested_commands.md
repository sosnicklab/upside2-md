
# Suggested Commands

## Installation & Build
```bash
# Install dependencies and compile C++ core (for Apple Silicon Macs)
./install_M1.sh

# Set up Python virtual environment
./install_python_env.sh

# Source environment variables
source source.sh
```

## Simulation Workflow
```bash
# Convert PDB to initial structure files
py/PDB_to_initial_structure.py input.pdb output_basename

# Create simulation configuration (.up file)
py/upside_config.py --output simulation.up --fasta output_basename.fasta ...

# Run constant temperature simulation
obj/upside --duration 1e7 --temperature 0.85 simulation.up

# Run replica exchange simulation
obj/upside --duration 1e7 --temperature 0.5,0.53,0.56 --swap-set 0-1 --swap-set 1-2 config_0.up config_1.up config_2.up
```

## Analysis
```bash
# List contents of .up file
py/attr_overview.py simulation.up

# Convert to VMD format
py/extract_vtf simulation.up simulation.vtf

# Load in Python with MDTraj
python -c "
import sys
sys.path.append('py')
import mdtraj_upside as mu
traj = mu.load_upside_traj('simulation.up')
print(f'Trajectory loaded with {len(traj)} frames')
"
```

## Utilities
```bash
# Check git status
git status

# List files in current directory
ls -la

# Find files with specific patterns
find . -name "*.cpp" -o -name "*.py"

# Search for patterns in files
grep -r "pattern" . --include="*.cpp" --include="*.py"
```
