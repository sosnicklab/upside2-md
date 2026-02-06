# Research Findings

## External References
- **MARTINI 3.0 Force Field**: https://www.martini.md/
- **CHARMM-GUI MARTINI Protocol**: https://charmm-gui.org/?doc=protocol/martini3

## API / Library Notes
- **Water Compressibility**: The compressibility of water is approximately 4.5e-5 bar⁻¹, which converts to ~2.1782 Å³/E_up. This is much smaller than the default value of 14.52118 Å³/E_up used for lipid systems.

## Barostat Configuration
- **Isotropic vs Semi-isotropic Coupling**:
  - Isotropic coupling (UPSIDE_NPT_SEMI=0): All axes scale uniformly in response to pressure changes
  - Semi-isotropic coupling (UPSIDE_NPT_SEMI=1): Lateral (x/y) and normal (z) axes scale independently
  - Water simulations should use isotropic coupling for uniform box behavior
- **Python Script Integration**: The workflow uses Python scripts (`prepare_martini.py` and `set_initial_position.py`) to prepare simulation inputs and update initial positions between stages.
- **HDF5 File Format**: The simulation configuration files are stored in HDF5 format, which allows for efficient storage and retrieval of large datasets.
- **Barostat Types**: The workflow uses Berendsen barostat for equilibration stages and Parrinello-Rahman barostat for production stage.
- **Softened Potentials**: Used during minimization and initial equilibration to avoid atomic clashes. Parameters include:
  - `lj_soften_alpha`: Controls Lennard-Jones potential softness (0=hard, 0.2=soft)
  - `slater_alpha`: Controls Coulomb potential softness (0=hard, 2.0=soft)
- **HDF5 Attribute Storage**: Softening parameters are stored as attributes in the HDF5 file under `/input/potential/martini_potential` and must be modified directly to change settings between stages.
