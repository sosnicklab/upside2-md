# GLY Ramachandran Phi-Symmetry: Literature Context and Fix

## The Problem

Glycine lacks a Cbeta atom, making its phi and -phi backbone conformations physically equivalent. Its Ramachandran distribution is therefore intrinsically symmetric about phi=0. However, PDB-derived statistical Ramachandran maps for GLY are asymmetric due to finite sampling noise, creating a spurious free-energy bias that favors one phi-quadrant over the other.

In Upside, this artifact caused TM helix 4 (TM4) of glpG to destabilize and exit the lipid bilayer during REMD. The spurious asymmetry in the GLY Rama map applied a net torque on TM4 backbone torsions that accumulated over the trajectory and drove the helix out. Symmetrizing the GLY map unconditionally in `write_rama_map_pot` fixed the problem.

## Literature Precedent

### Rosetta

Rosetta independently identified and documented this exact problem. They introduced the flag `-symmetric_gly_tables` (around 2016) because their PDB-derived GLY `rama` and `rama_prepro` tables are asymmetric from sampling noise, creating a bias toward D-amino acid conformations. Their docs state:

> "By default, the gly table is asymmetric because it is based on statistics from the PDB, which disproportionately puts glycine in the D-amino acid region of Ramachandran space."

The flag is off by default, meaning most Rosetta users run with the spuriously asymmetric table unless they opt in. Sources:
- Rosetta Ramachandran Class Reference
- Rosetta `simple_cycpep_predict` documentation

### Atomistic Force Fields (CHARMM36, AMBER ff14SB/ff19SB)

These sidestep the issue by deriving backbone torsion surfaces from quantum chemistry on the glycine dipeptide rather than PDB statistics. The QM calculation naturally produces a symmetric surface, so no explicit symmetrization is needed. The problem only arises when a statistical (PDB-derived) potential is used.

### Coarse-Grained Models (MARTINI, UNRES, CABS, PRIMO)

None of these explicitly discuss GLY phi-symmetry in their published parameterizations.

### TM Context

No prior paper was found attributing transmembrane helix instability or bilayer exit to a GLY Ramachandran asymmetry artifact. The glpG TM4 case appears to be the first documented instance of this mechanism in a membrane simulation context.

## Fix in Upside

Symmetrize ALL GLY Ramachandran maps unconditionally in `write_rama_map_pot`. Any phi-range condition or `alphaR > alphaL` guard silently misses GLY residues and breaks TM helix stability. See memory entry `feedback_gly_rama_symmetrization.md` for the implementation rule.

## Key References

- Shelar & Chattopadhyay (2005). *The Ramachandran plots of glycine and pre-proline.* BMC Struct Biol 5:14. Full analysis of GLY's distinct symmetric map. https://pmc.ncbi.nlm.nih.gov/articles/PMC1201153/
- Baruch-Shpigler et al. (2017). *Chiral Ramachandran Plots I: Glycine.* Biochemistry. Confirms GLY occupies both phi-quadrants essentially equally in high-resolution structures. https://pubs.acs.org/doi/10.1021/acs.biochem.7b00525
- Lovell et al. (2003). *Structure validation by Calpha geometry: phi, psi and Cbeta deviation.* Proteins. Standard reference establishing GLY's symmetric Ramachandran as the validation baseline.
- Hovmoller, Zhou, Ohlson (2002). *Conformations of amino acids in proteins.* Acta Crystallogr D. GLY's two-fold symmetric distribution documented.
