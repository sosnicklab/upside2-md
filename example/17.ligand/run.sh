#!/bin/bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")"/../.. && pwd)"
EXAMPLE_ROOT="$PROJECT_ROOT/example/17.ligand"

source "$PROJECT_ROOT/.venv/bin/activate"
source "$PROJECT_ROOT/source.sh"

cd "$EXAMPLE_ROOT"
python3 0.prepare.py
python3 1.run.py
python3 2.validate.py
python3 3.visualize.py
python3 4.compare_forcefield.py
