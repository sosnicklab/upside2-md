#!/usr/bin/env bash

if [ -z "$UPSIDE_HOME" ]; then
    echo "Error: UPSIDE_HOME is not set."
    exit 1
fi

PYTHON="$UPSIDE_HOME/.venv/bin/python3"
SCRIPT="$UPSIDE_HOME/tests/validate_method.py"

#interaction_graph.h
"$PYTHON" "$SCRIPT" --method=change_cache_buffer --args=args.new_buffer --json="$UPSIDE_HOME/tests/records/line_182_capture.json"