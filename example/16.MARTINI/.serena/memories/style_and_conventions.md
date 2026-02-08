# Style and Conventions for 16.MARTINI Project

## Programming Languages
- **Shell Scripts**: Bash (shebang: #!/bin/bash)
- **Python Scripts**: Python 3.x

## Shell Script Style
- Use `.sh` extension for shell scripts
- Shebang line required: `#!/bin/bash`
- Use lowercase filenames with underscores
- Comment sections for stage documentation
- Error handling with `set -euo pipefail`
- Explicit parameterization with environment variables

## Python Script Style
- Use `.py` extension for Python scripts
- Follow Python PEP8 style guidelines
- Use lowercase filenames with underscores
- Docstrings for functions and modules
- Error handling with try-except blocks
- Command-line argument parsing with `argparse`

## Project Structure
- `ff3.00/`: MARTINI force field files
- `pdb/`: PDB structure files
- `*.sh`: Shell scripts for simulation execution
- `*.py`: Python scripts for preparation and analysis
- `task_plan.md`: Project goal and execution plan
- `findings.md`: Technical research findings
- `progress.md`: Execution log

## Configuration Management
- Simulation parameters passed via command-line arguments or environment variables
- HDF5 files used for storing simulation configurations
- Shell scripts use explicit parameterization for reproducibility
