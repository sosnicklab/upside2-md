#!.venv/bin/python
import mdtraj as md

# Step 1: Load the full structure from PDB
full_traj = md.load('charmm-gui-4431828778/gromacs/step5_charmm2gmx.pdb')
print(f"Full PDB: {full_traj.n_atoms} atoms, {full_traj.n_residues} residues")

"""

# Step 2: Load individual molecule topologies from .gro files
molecule_topos = {}
gro_files = ['mol1.gro', 'mol2.gro', 'mol3.gro']  # or use glob.glob('*.gro')

for gro_file in gro_files:
    mol_traj = md.load(gro_file)
    molecule_topos[gro_file] = mol_traj.topology
    print(f"{gro_file}: {mol_traj.n_atoms} atoms")

# (Optional) Step 3: Match or process topology info
# You could loop through residues in full_traj.topology and annotate/match with molecule_topos

# Example: list atoms in first molecule topology
topo = molecule_topos['mol1.gro']
for atom in topo.atoms:
    print(f"{atom.index}: {atom.name} in {atom.residue}")
"""