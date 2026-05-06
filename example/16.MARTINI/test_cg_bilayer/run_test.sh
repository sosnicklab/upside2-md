#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

cd "$REPO_ROOT"
source .venv/bin/activate
source source.sh

set -euo pipefail

python3 example/16.MARTINI/test_cg_bilayer/run_test.py "$@"
