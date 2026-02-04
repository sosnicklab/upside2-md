# Research Findings

## External References
- **MARTINI 3.0 Force Field**: https://www.martini.md/
- **CHARMM-GUI MARTINI Protocol**: https://charmm-gui.org/?doc=protocol/martini3

## API / Library Notes
- **Python Script Integration**: The workflow uses Python scripts (`prepare_martini.py` and `set_initial_position.py`) to prepare simulation inputs and update initial positions between stages.
- **HDF5 File Format**: The simulation configuration files are stored in HDF5 format, which allows for efficient storage and retrieval of large datasets.
- **Barostat Types**: The workflow uses Berendsen barostat for equilibration stages and Parrinello-Rahman barostat for production stage.
- **Softened Potentials**: Used during minimization and initial equilibration to avoid atomic clashes. Parameters include:
  - `lj_soften_alpha`: Controls Lennard-Jones potential softness (0=hard, 0.2=soft)
  - `slater_alpha`: Controls Coulomb potential softness (0=hard, 2.0=soft)
- **HDF5 Attribute Storage**: Softening parameters are stored as attributes in the HDF5 file under `/input/potential/martini_potential` and must be modified directly to change settings between stages.
