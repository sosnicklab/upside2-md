#!/usr/bin/env python3
"""Rewrite every replica's Energy.npy with the full coupled-system potential, referenced to its mean.

Two things are fixed here at once, both required before helpers/calc_hdx_ht.py can reweight a hybrid
dry-MARTINI REMD run.

The energy: example/00.AnalysisScripts/README.md fixes the contract -- for a hybrid trajectory
Energy.npy "must be the full coupled-system potential, not the protein-only analysis energy".
martini_hdx_project.py re-scores /output/potential to the protein-only energy, so
get_info_from_upside_traj.py reads that instead. The coupled potential is taken from the trimmed
trajectories, whose frames are the ones the projection was built from, so the arrays stay aligned
frame for frame.

The reference: the coupled potential of a protein plus a DDM micelle is around -7.6e3 E_up, giving
reduced potentials near -1.2e4, and exp(-u) at that magnitude overflows. Measured on
glpG-RKRK-79HIS_S115T, MBAR then never leaves its initial guess (f_k spread 0.000, neighbour overlap
0.0000) and the reweighting weights underflow to zero at every target temperature but the lowest rung.
Subtracting one constant C from every replica's energy is exact, not an approximation: it shifts f_k by
-beta_k*C and multiplies both the numerator and the denominator of every reweighting weight by
exp(beta_target*C), which cancels in the normalisation. With C set to the pooled mean the same run
solves (f_k spread 71.07, neighbour overlap 0.115-0.128, all frames carrying weight at every target).
C must be one constant shared by all replicas, so all of them are written in a single pass.
"""
import argparse
import os

import h5py
import numpy as np


def hybrid_potential(path):
    with h5py.File(path, "r") as h5:
        potential = np.asarray(h5["/output/potential"][:]).reshape(-1)
    if not np.isfinite(potential).all():
        raise SystemExit("%s has non-finite potentials" % path)
    return potential


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--trimmed-dir", required=True, help="directory of trimmed hybrid trajectories")
    parser.add_argument("--results-dir", required=True, help="directory holding the per-replica arrays")
    parser.add_argument("--pdb-id", required=True, help="variant name")
    parser.add_argument("--sim-id", required=True, help="simulation tag used in the array filenames")
    parser.add_argument("--n-replica", type=int, required=True, help="number of replicas")
    args = parser.parse_args()

    trimmed = ["%s/%s.run.%d.up" % (args.trimmed_dir, args.pdb_id, i) for i in range(args.n_replica)]
    targets = ["%s/%s_%s_%d_Energy.npy" % (args.results_dir, args.pdb_id, args.sim_id, i)
               for i in range(args.n_replica)]
    for path in trimmed + targets:
        if not os.path.isfile(path):
            raise SystemExit("missing %s" % path)

    potentials = [hybrid_potential(path) for path in trimmed]
    for i, (potential, target) in enumerate(zip(potentials, targets)):
        existing = np.load(target).reshape(-1).size
        if potential.size != existing:
            raise SystemExit("replica %d: hybrid potential has %d frames but Energy.npy has %d"
                             % (i, potential.size, existing))

    reference = float(np.concatenate(potentials).mean())
    for potential, target in zip(potentials, targets):
        np.save(target, potential - reference)
    print("Energy.npy <- coupled-system potential, reference %.4f E_up subtracted from %d replicas"
          % (reference, args.n_replica), flush=True)


if __name__ == "__main__":
    main()
