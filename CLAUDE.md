# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview
Upside is a molecular dynamics simulation package for protein folding and conformational dynamics. It combines a fast C++ core with Python scripts for configuration and analysis.

## Common Commands

### Environment Setup
Crucial: You must run these commands from the project root before running anything in this project:
```bash
source .venv/bin/activate
source source.sh

```

### Installation & Build
```bash
# Install dependencies and compile C++ core
./install_M1.sh
./install_python_env.sh
```

### Upside Unit Conversions

| Quantity | Upside Unit | Standard Equivalent |
| :--- | :--- | :--- |
| **Energy** | 1 E_up | 2.914952774272 kJ/mol |
| **Length** | 1 Angstrom | 1 Angstrom |
| **Mass** | 1 m_up | 12 g/mol |
| **Temperature** | 1.0 T_up | 350.588235 Kelvin |
| **Pressure** | 0.000020933215 E_up / (Angstrom^3) | 1 atm |
| **Pressure** | 0.000020659477 E_up /(Angstrom^3) | 1 bar |
| **Compressibility** | 14.521180763676 Angstrom^3 / E_up | 3e-4 bar^(-1) |
