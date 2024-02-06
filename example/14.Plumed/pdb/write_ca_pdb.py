#!/usr/bin/env python
import MDAnalysis as mda

univ = mda.Universe("chig.pdb")
ca = univ.select_atoms("name CA")
with mda.Writer("chig.ca.pdb", ca.n_atoms) as f:
    f.write(ca)
