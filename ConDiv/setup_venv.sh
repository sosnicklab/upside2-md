#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
VENV_DIR="${SCRIPT_DIR}/venv"

if [ ! -d "$VENV_DIR" ]; then
  python3 -m venv "$VENV_DIR"
fi

source "$VENV_DIR/bin/activate"
python3 -m pip install --upgrade pip setuptools wheel
python3 -m pip install numpy scipy tables mdtraj

echo "ConDiv local venv ready: $VENV_DIR"
