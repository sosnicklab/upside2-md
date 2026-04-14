#!/bin/bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

REQUIRED_PYTHON_MM="3.11"
VENV_DIR="${VENV_DIR:-$ROOT_DIR/.venv}"
HOST_OS="$(uname -s)"

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

resolve_executable() {
    local candidate="$1"

    if [[ "$candidate" == */* ]]; then
        [ -x "$candidate" ] || return 1
        printf '%s\n' "$candidate"
        return 0
    fi

    command -v "$candidate" 2>/dev/null
}

python_major_minor() {
    local python_bin="$1"
    "$python_bin" -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")' 2>/dev/null
}

discover_python_bin() {
    local candidates=()
    local candidate resolved candidate_mm

    case "$HOST_OS" in
        Darwin)
            candidates=(
                python3.11
                /opt/homebrew/bin/python3.11
                /usr/local/bin/python3.11
                python3
            )
            ;;
        Linux)
            candidates=(
                python3.11
                /usr/bin/python3.11
                /usr/local/bin/python3.11
                python3
            )
            ;;
        *)
            candidates=(python3.11 python3)
            ;;
    esac

    for candidate in "${candidates[@]}"; do
        resolved="$(resolve_executable "$candidate")" || continue
        candidate_mm="$(python_major_minor "$resolved" || true)"
        if [ "$candidate_mm" = "$REQUIRED_PYTHON_MM" ]; then
            printf '%s\n' "$resolved"
            return 0
        fi
    done

    return 1
}

prepend_space_var() {
    local var_name="$1"
    local new_value="$2"
    local current_value="${!var_name:-}"

    if [ -z "$current_value" ]; then
        export "$var_name=$new_value"
    else
        export "$var_name=$new_value $current_value"
    fi
}

prepend_path_var() {
    local var_name="$1"
    local new_value="$2"
    local current_value="${!var_name:-}"

    if [ -z "$current_value" ]; then
        export "$var_name=$new_value"
    else
        export "$var_name=$new_value:$current_value"
    fi
}

add_prefix_build_hints() {
    local prefix="$1"

    [ -d "$prefix/include" ] && prepend_space_var CPPFLAGS "-I$prefix/include"
    [ -d "$prefix/lib" ] && prepend_space_var LDFLAGS "-L$prefix/lib"
    [ -d "$prefix/lib/pkgconfig" ] && prepend_path_var PKG_CONFIG_PATH "$prefix/lib/pkgconfig"
    return 0
}

configure_macos_build_env() {
    local prefix

    if [ "$HOST_OS" != "Darwin" ]; then
        return
    fi

    if ! command -v brew >/dev/null 2>&1; then
        echo "WARNING: Homebrew not found. Native-extension builds may fail on macOS." >&2
        return
    fi

    if prefix="$(brew --prefix hdf5 2>/dev/null)"; then
        export HDF5_DIR="${HDF5_DIR:-$prefix}"
        export HDF5_ROOT="${HDF5_ROOT:-$prefix}"
        add_prefix_build_hints "$prefix"
    fi

    if prefix="$(brew --prefix c-blosc2 2>/dev/null)"; then
        export BLOSC2_DIR="${BLOSC2_DIR:-$prefix}"
        add_prefix_build_hints "$prefix"
    fi

    if prefix="$(brew --prefix lzo 2>/dev/null)"; then
        export LZO_DIR="${LZO_DIR:-$prefix}"
        add_prefix_build_hints "$prefix"
    fi

    if prefix="$(brew --prefix bzip2 2>/dev/null)"; then
        export BZIP2_DIR="${BZIP2_DIR:-$prefix}"
        add_prefix_build_hints "$prefix"
    fi

    return 0
}

if [ -n "${PYTHON_BIN:-}" ]; then
    if ! PYTHON_BIN="$(resolve_executable "$PYTHON_BIN")"; then
        echo "ERROR: Requested Python interpreter not found: ${PYTHON_BIN}" >&2
        exit 1
    fi
else
    if ! PYTHON_BIN="$(discover_python_bin)"; then
        echo "ERROR: Could not find a Python $REQUIRED_PYTHON_MM interpreter. Install python3.11 or set PYTHON_BIN=/path/to/python3.11" >&2
        exit 1
    fi
fi

echo "Setting up Python environment for UPSIDE2 in $VENV_DIR"
echo "Host OS: $HOST_OS"
echo "Using Python interpreter: $PYTHON_BIN"

host_python_mm="$(python_major_minor "$PYTHON_BIN" || true)"
if [ "$host_python_mm" != "$REQUIRED_PYTHON_MM" ]; then
    echo "ERROR: install_python_env.sh only supports Python $REQUIRED_PYTHON_MM. $PYTHON_BIN resolves to Python ${host_python_mm:-unknown}" >&2
    exit 1
fi

configure_macos_build_env

recreate_venv=false
if [ -d "$VENV_DIR" ]; then
    if [ ! -x "$VENV_DIR/bin/python" ]; then
        recreate_venv=true
    else
        current_venv_mm="$(python_major_minor "$VENV_DIR/bin/python" || true)"
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

if [ -x "$VENV_DIR/bin/python" ]; then
    VENV_PYTHON="$VENV_DIR/bin/python"
elif [ -x "$VENV_DIR/bin/python3" ]; then
    VENV_PYTHON="$VENV_DIR/bin/python3"
else
    echo "ERROR: No Python interpreter found inside $VENV_DIR/bin" >&2
    exit 1
fi

"$VENV_PYTHON" -m pip install --upgrade --prefer-binary pip setuptools wheel
"$VENV_PYTHON" -m pip install --upgrade --prefer-binary "${CORE_PACKAGES[@]}"

failed_optional_packages=()
for package in "${OPTIONAL_PACKAGES[@]}"; do
    if ! "$VENV_PYTHON" -m pip install --upgrade --prefer-binary "$package"; then
        failed_optional_packages+=("$package")
        echo "WARNING: Optional package failed to install: $package" >&2
    fi
done

VERIFY_CACHE_ROOT="${TMPDIR:-/tmp}/upside-python-install-$RANDOM-$$"
mkdir -p "$VERIFY_CACHE_ROOT/mpl" "$VERIFY_CACHE_ROOT/cache"
trap 'rm -rf "$VERIFY_CACHE_ROOT"' EXIT

MPLCONFIGDIR="$VERIFY_CACHE_ROOT/mpl" XDG_CACHE_HOME="$VERIFY_CACHE_ROOT/cache" "$VENV_PYTHON" - <<'PY'
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

if [ "${#failed_optional_packages[@]}" -gt 0 ]; then
    echo "Optional packages not installed: ${failed_optional_packages[*]}"
fi

echo "Python environment setup completed."
echo "Virtual environment available at: $VENV_DIR"
