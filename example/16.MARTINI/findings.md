# Research Findings

## External References
- **MARTINI 3.0 Force Field**: https://www.martini.md/
- **CHARMM-GUI MARTINI Protocol**: https://charmm-gui.org/?doc=protocol/martini3

## API / Library Notes
- **Python Script Integration**: The workflow uses Python scripts (`prepare_martini.py` and `set_initial_position.py`) to prepare simulation inputs and update initial positions between stages.
- **HDF5 File Format**: The simulation configuration files are stored in HDF5 format, which allows for efficient storage and retrieval of large datasets.
- **Barostat Types**: The workflow uses Berendsen barostat for equilibration stages (6.2-6.6) and Parrinello-Rahman barostat for production stage (7.0).
