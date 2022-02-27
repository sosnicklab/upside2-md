#!/usr/bin/env python
import sys, os
import mdtraj as md
import mdtraj_upside as mu

ref_file = sys.argv[1]
re = mu.load_upside_ref(ref_file) # load the structure of the input

traj_file = sys.argv[2]
tr = mu.load_upside_traj(traj_file) 

sele = tr.top.select("name CA")
rmsd = 10.*md.rmsd(tr, re, atom_indices=sele)

for r in rmsd:
    print(r)
