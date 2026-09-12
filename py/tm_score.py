#!/usr/bin/env python
"""TM-score against a reference structure, for a fixed residue correspondence.

TM = max over superpositions of (1/L_norm) * sum_i 1/(1 + (d_i/d0)^2), with
d0 = 1.24*(L_norm-15)^(1/3) - 1.8 floored at 0.5 (Zhang & Skolnick, Proteins 57:702, 2004).

This is TM-score, not TM-align: the residue correspondence is fixed at i <-> i and only the
superposition is optimised. That is the right quantity when model and reference are the same
protein, which is the case for a folding benchmark scored against its own native structure.
TM-align is a different number because it is also free to shift the alignment, and it reads
*higher* on poorly folded structures for that reason; do not use it as a reference below
TM ~ 0.5.

The maximum is located as in the reference program: seed the superposition from contiguous
fragments of decreasing length at sliding offsets, then iteratively re-superimpose on the
residues inside a distance cutoff until the selected set stops changing, keeping the best
score over all seeds. Several starting cutoffs are tried because that iteration is only
locally optimal.

Validated 2026-09-12 against tmtools 0.3.0 (a binding to the reference TM-align C++ code):

* d0 reproduces the published formula exactly (L=76 -> 3.081).
* A structure against itself gives 1.000000, and the score is invariant to rigid rotation
  and translation (1.000000).
* Over 19 synthetic pairs with TM > 0.5 the largest disagreement with the reference is
  0.007 and the mean is 0.0006.
* The seed density does not change the answer in the folded regime: step = 16, 8, 4 and 1
  agree to <= 1e-4 for TM > 0.5, so the optimiser is converged rather than merely stopping.
  Below TM ~ 0.2 the step does matter at the ~0.01 level; use step = 1 if the unfolded tail
  of a distribution has to be quantitative.

`step` trades seed density for speed: step = 2 costs about 65 ms per 76-residue structure
and is what the ff3.0 benchmark used.
"""
import numpy as np


def d0_of(n_residue):
    """Length scale of the TM-score sum, floored as the reference program floors it."""
    return max(0.5, 1.24 * (max(n_residue, 16) - 15) ** (1.0 / 3.0) - 1.8)


def _superpose(mobile, target):
    """Rotation and translation carrying mobile onto target, both (n, 3)."""
    mobile_center, target_center = mobile.mean(0), target.mean(0)
    u, _, vt = np.linalg.svd((mobile - mobile_center).T @ (target - target_center))
    chirality = np.sign(np.linalg.det(vt.T @ u.T))
    rotation = vt.T @ np.diag([1.0, 1.0, chirality]) @ u.T
    return rotation, target_center - rotation @ mobile_center


def tm_score(model, reference, length_norm=None, step=2):
    """TM-score of model against reference, both (L, 3) CA coordinates in the same order.

    length_norm defaults to L, which is the convention when the reference is the native
    structure of the same protein.
    """
    model = np.ascontiguousarray(model, dtype=float)
    reference = np.ascontiguousarray(reference, dtype=float)
    n_residue = len(reference)
    if len(model) != n_residue:
        raise ValueError('model has %d residues, reference %d; TM-score needs a 1:1 '
                         'correspondence' % (len(model), n_residue))
    d0 = d0_of(length_norm or n_residue)
    d0_sq = d0 * d0

    seed_lengths = []
    length = n_residue
    while length >= 4:
        seed_lengths.append(length)
        length //= 2
    if 4 not in seed_lengths:
        seed_lengths.append(4)

    best = 0.0
    for seed_length in seed_lengths:
        offset = max(1, step) if seed_length < n_residue else 1
        for start in range(0, n_residue - seed_length + 1, offset):
            seed = np.arange(start, start + seed_length)
            for cutoff0 in (d0 - 1.0, d0, d0 + 1.0, d0 + 2.0):
                selected, previous = seed, None
                for _ in range(30):
                    rotation, translation = _superpose(model[selected], reference[selected])
                    delta = (model @ rotation.T + translation) - reference
                    d_sq = np.einsum('ij,ij->i', delta, delta)
                    score = np.mean(1.0 / (1.0 + d_sq / d0_sq))
                    if score > best:
                        best = score
                    # grow the cutoff until at least four residues survive, else the
                    # superposition is undetermined and this seed is finished
                    cutoff = max(cutoff0, 0.5)
                    while True:
                        keep = np.where(d_sq < cutoff * cutoff)[0]
                        if len(keep) >= 4 or cutoff > 25.0:
                            break
                        cutoff += 0.5
                    if len(keep) < 4:
                        break
                    if previous is not None and np.array_equal(keep, previous):
                        break
                    previous, selected = selected, keep
    return best


def main():
    import sys
    import mdtraj as md
    import mdtraj_upside as mu

    if len(sys.argv) < 3:
        sys.exit('usage: tm_score.py <reference.up> <trajectory.up>\n'
                 '  prints the TM-score of every frame against the reference')
    reference = mu.load_upside_ref(sys.argv[1])
    traj = mu.load_upside_traj(sys.argv[2])
    ca = traj.top.select('name CA')
    ref_ca = 10.0 * reference.xyz[0][reference.top.select('name CA')]
    for frame in traj.xyz:
        print(tm_score(10.0 * frame[ca], ref_ca))


if __name__ == '__main__':
    main()
