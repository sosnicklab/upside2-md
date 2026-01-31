# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview
Upside is a molecular dynamics simulation package for protein folding and conformational dynamics. It combines a fast C++ core with Python scripts for configuration and analysis.

## Common Commands

### Installation & Build
```bash
# Install dependencies and compile C++ core
./install_M1.sh
./install_python_env.sh
```

### Environment Setup
Crucial: You must run these commands from the project root before running any scripts or simulations:
```bash
source .venv/bin/activate
source source.sh

```