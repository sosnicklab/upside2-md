#!/bin/bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

PYTHON_BIN="${PYTHON_BIN:-python3.11}"
REQUIRED_PYTHON_MM="3.11"
VENV_DIR="$ROOT_DIR/.venv"

CORE_PACKAGES=(
    h5py
    tables
    matplotlib
    mdtraj
    pymbar
    pandas
    ProDy
    scikit-learn
    jax
)

OPTIONAL_PACKAGES=(
    colorcet
    pyhdx==0.4.3
    "hdxms-datasets<0.2"
)

echo "Setting up Python environment for UPSIDE2 in $VENV_DIR"

if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
    echo "ERROR: Requested Python interpreter not found: $PYTHON_BIN" >&2
    exit 1
fi

host_python_mm="$("$PYTHON_BIN" -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')"
if [ "$host_python_mm" != "$REQUIRED_PYTHON_MM" ]; then
    echo "ERROR: install_python_env.sh only supports Python $REQUIRED_PYTHON_MM. $PYTHON_BIN resolves to Python $host_python_mm" >&2
    exit 1
fi

recreate_venv=false
if [ -d "$VENV_DIR" ]; then
    if [ ! -x "$VENV_DIR/bin/python" ]; then
        recreate_venv=true
    else
        current_venv_mm="$("$VENV_DIR/bin/python" -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')"
        if [ "$current_venv_mm" != "$REQUIRED_PYTHON_MM" ]; then
            recreate_venv=true
            echo "Existing .venv uses Python $current_venv_mm; recreating it with $PYTHON_BIN"
        fi
    fi
fi

if [ "$recreate_venv" = true ]; then
    rm -rf "$VENV_DIR"
fi

if [ ! -d "$VENV_DIR" ]; then
    "$PYTHON_BIN" -m venv "$VENV_DIR"
fi

source "$VENV_DIR/bin/activate"

python -m pip install --upgrade pip setuptools wheel
python -m pip install --upgrade "${CORE_PACKAGES[@]}" "${OPTIONAL_PACKAGES[@]}"

python - <<'PY'
modules = [
    "h5py",
    "tables",
    "matplotlib",
    "mdtraj",
    "pymbar",
    "pandas",
    "prody",
    "sklearn",
    "jax",
]
for module in modules:
    __import__(module)
print("Core scientific stack import check passed.")
print(f"Python executable: {__import__('sys').executable}")
print(f"Python version: {__import__('sys').version.split()[0]}")

try:
    import pyhdx  # noqa: F401
except Exception as exc:
    print(f"WARNING: PyHDX import check failed: {exc}")
else:
    print("PyHDX import check passed.")
PY

deactivate

echo "Python environment setup completed."
echo "Virtual environment available at: $VENV_DIR"
