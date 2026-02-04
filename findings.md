# Research Findings

## External References
- CHARMM-GUI MARTINI Protocol:
  - Steps 6.0-7.0 for lipid bilayer equilibration
  - Requires Parrinello-Rahman barostat for production runs
  - Uses position restraints on lipid headgroups during equilibration

- Parrinello-Rahman Barostat:
  - Extended Lagrangian method for NPT ensemble
  - Box matrix h evolves according to: dh/dt = V/W * (P - P_target) * h
  - W is barostat mass parameter (related to tau_p)
  - Provides correct volume fluctuations unlike Berendsen
  - Semi-isotropic: h_xx = h_yy (lateral), h_zz (normal) independent
  - Standard implementation in GROMACS and LAMMPS

## API / Library Notes
- UPSIDE Barostat:
  - Currently only supports Berendsen algorithm (lines 229-261 in src/box.cpp)
  - Located in src/box.cpp
  - Needs Parrinello-Rahman implementation for proper ensemble properties
  - Current structure: BarostatSettings and BarostatState in box.h

- UPSIDE MARTINI:
  - Located in src/martini.cpp
  - Needs position restraint potential node for headgroup restraints
