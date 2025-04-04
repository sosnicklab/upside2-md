# align_pdbs.py

import argparse
import MDAnalysis as mda
from MDAnalysis.analysis import align

def main(all_atom_pdb, martini_pdb, output_file, update_target):
    # Load the structures
    aa = mda.Universe(all_atom_pdb)
    cg = mda.Universe(martini_pdb)

    # Select backbone or representative atoms (simplified matching)
    # Modify these selections based on your mapping
    aa_sel = aa.select_atoms("name CA or name BB or backbone")
    cg_sel = cg.select_atoms("name BB")

    if len(aa_sel) != len(cg_sel):
        raise ValueError(f"Selections must be same size. Got {len(aa_sel)} (AA) and {len(cg_sel)} (CG).")

    if update_target == "martini":
        # Align MARTINI to all-atom
        aligner = align.AlignTraj(cg, aa, select="name BB", in_memory=True)
        aligner.run()
        cg.atoms.write(output_file)
    elif update_target == "all_atom":
        # Align all-atom to MARTINI
        aligner = align.AlignTraj(aa, cg, select="name CA or name BB or backbone", in_memory=True)
        aligner.run()
        aa.atoms.write(output_file)
    else:
        raise ValueError("update_target must be 'martini' or 'all_atom'.")

    print(f"Alignment complete. Updated structure written to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="RMSD align MARTINI and all-atom PDBs.")
    parser.add_argument("--all_atom_pdb", required=True, help="Path to all-atom PDB file")
    parser.add_argument("--martini_pdb", required=True, help="Path to MARTINI coarse-grained PDB file")
    parser.add_argument("--output", required=True, help="Path to write aligned PDB")
    parser.add_argument("--update_target", choices=["martini", "all_atom"], required=True,
                        help="Which structure to update coordinates for")

    args = parser.parse_args()
    main(args.all_atom_pdb, args.martini_pdb, args.output, args.update_target)
