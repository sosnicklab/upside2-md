#!/usr/bin/env python3

import math
import os
import sys

import numpy as np
import tables as tb

import charmm_import as ci


SIM_ID = "smoke"
SCALAR_ABS_TOL = 1e-3
SCALAR_REL_TOL = 1e-5
ARRAY_ABS_TOL = 1e-3
ARRAY_REL_TOL = 2e-5


def example_root():
    return os.path.abspath(os.path.dirname(__file__))


def repo_root():
    return os.path.abspath(os.path.join(example_root(), "..", ".."))


def upside_home():
    return os.environ.get("UPSIDE_HOME", repo_root())


def ensure_loader_compat():
    obj_dir = os.path.join(upside_home(), "obj")
    dylib = os.path.join(obj_dir, "libupside.dylib")
    so = os.path.join(obj_dir, "libupside.so")
    if os.path.exists(dylib) and not os.path.exists(so):
        os.symlink(dylib, so)


def load_upside_engine():
    ensure_loader_compat()
    py_path = os.path.join(upside_home(), "py")
    if py_path not in sys.path:
        sys.path.insert(0, py_path)
    import upside_engine as ue

    return ue


def run_path():
    return os.path.join(example_root(), "outputs", SIM_ID, "4w52_ligand.run.up")


def latest_output_group(handle):
    idx = 0
    latest = None
    while "output_previous_%d" % idx in handle.root:
        latest = handle.get_node("/output_previous_%d" % idx)
        idx += 1
    if "output" in handle.root:
        latest = handle.get_node("/output")
    return latest


def try_openmm_status():
    try:
        import openmm  # noqa: F401

        return "available"
    except Exception as exc:
        return "unavailable (%s)" % type(exc).__name__


def quat_to_rot(quat):
    a, b, c, d = quat
    return np.array([
        [a * a + b * b - c * c - d * d, 2.0 * b * c - 2.0 * a * d, 2.0 * b * d + 2.0 * a * c],
        [2.0 * b * c + 2.0 * a * d, a * a - b * b + c * c - d * d, 2.0 * c * d - 2.0 * a * b],
        [2.0 * b * d - 2.0 * a * c, 2.0 * c * d + 2.0 * a * b, a * a - b * b - c * c + d * d],
    ], dtype="f8")


def apply_affine(affine, local_pos):
    rot = quat_to_rot(affine[3:])
    return local_pos.dot(rot.T) + affine[:3][None, :]


def charmm_pair_energy(r2, rmin, epsilon, qq, cutoff):
    if r2 >= cutoff * cutoff:
        return 0.0
    safe_r2 = max(float(r2), 1e-6)
    inv_r2 = 1.0 / safe_r2
    rmin2_over_r2 = (float(rmin) * float(rmin)) * inv_r2
    sr6 = rmin2_over_r2 * rmin2_over_r2 * rmin2_over_r2
    sr12 = sr6 * sr6
    lj = float(epsilon) * (sr12 - 2.0 * sr6)
    coulomb = float(qq) * math.sqrt(inv_r2)
    return lj + coulomb


def dihedral_angle(r1, r2, r3, r4):
    f = r1 - r2
    g = r2 - r3
    h = r4 - r3
    a = np.cross(f, g)
    b = np.cross(h, g)
    c = np.cross(b, a)
    gmag = np.linalg.norm(g)
    return math.atan2(float(np.dot(c, g)), float(np.dot(a, b) * gmag))


def wrap_angle(angle):
    while angle > math.pi:
        angle -= 2.0 * math.pi
    while angle < -math.pi:
        angle += 2.0 * math.pi
    return angle


def compare_scalar(name, ref_value, engine_value):
    diff = abs(ref_value - engine_value)
    ok = np.isclose(ref_value, engine_value, atol=SCALAR_ABS_TOL, rtol=SCALAR_REL_TOL)
    print("%-32s ref=%12.6f engine=%12.6f diff=%10.3e" % (name, ref_value, engine_value, diff))
    if not ok:
        raise RuntimeError("%s mismatch: ref=%r engine=%r" % (name, ref_value, engine_value))


def compare_array(name, ref_value, engine_value):
    ref_value = np.asarray(ref_value, dtype="f8")
    engine_value = np.asarray(engine_value, dtype="f8")
    diff = np.abs(ref_value - engine_value)
    scale = np.maximum.reduce([np.ones_like(diff), np.abs(ref_value), np.abs(engine_value)])
    rel = diff / scale
    max_abs = float(diff.max()) if diff.size else 0.0
    max_rel = float(rel.max()) if rel.size else 0.0
    print("%-32s shape=%-12s max_abs=%10.3e max_rel=%10.3e" % (name, str(ref_value.shape), max_abs, max_rel))
    if not np.allclose(ref_value, engine_value, atol=ARRAY_ABS_TOL, rtol=ARRAY_REL_TOL):
        idx = np.unravel_index(np.argmax(rel), rel.shape)
        raise RuntimeError(
            "%s mismatch at %s: ref=%r engine=%r abs=%r rel=%r"
            % (name, idx, ref_value[idx], engine_value[idx], diff[idx], rel[idx])
        )


def scalar_output(engine, node_name):
    return float(engine.get_output(node_name).reshape(-1)[0])


def bond_energy(group, ligand_pos):
    energy = 0.0
    atom1 = group.atom1[:]
    atom2 = group.atom2[:]
    length = group.length[:]
    spring = group.spring[:]
    for a, b, l0, k in zip(atom1, atom2, length, spring):
        disp = ligand_pos[a] - ligand_pos[b]
        delta = np.linalg.norm(disp) - float(l0)
        energy += 0.5 * float(k) * delta * delta
    return energy


def angle_energy(group, ligand_pos):
    energy = 0.0
    for a, b, c, theta0, spring in zip(
        group.atom1[:],
        group.atom2[:],
        group.atom3[:],
        group.theta0[:],
        group.spring[:],
    ):
        x1 = ligand_pos[a] - ligand_pos[b]
        x2 = ligand_pos[c] - ligand_pos[b]
        x1h = x1 / max(np.linalg.norm(x1), 1e-6)
        x2h = x2 / max(np.linalg.norm(x2), 1e-6)
        theta = math.acos(float(np.clip(np.dot(x1h, x2h), -1.0, 1.0)))
        delta = theta - float(theta0)
        energy += 0.5 * float(spring) * delta * delta
    return energy


def torsion_energy(group, ligand_pos):
    energy = 0.0
    for a, b, c, d, periodicity, phase, amplitude in zip(
        group.atom1[:],
        group.atom2[:],
        group.atom3[:],
        group.atom4[:],
        group.periodicity[:],
        group.phase[:],
        group.amplitude[:],
    ):
        phi = dihedral_angle(ligand_pos[a], ligand_pos[b], ligand_pos[c], ligand_pos[d])
        energy += float(amplitude) * (1.0 + math.cos(float(periodicity) * phi - float(phase)))
    return energy


def improper_energy(group, ligand_pos):
    energy = 0.0
    for a, b, c, d, theta0, spring in zip(
        group.atom1[:],
        group.atom2[:],
        group.atom3[:],
        group.atom4[:],
        group.theta0[:],
        group.spring[:],
    ):
        phi = dihedral_angle(ligand_pos[a], ligand_pos[b], ligand_pos[c], ligand_pos[d])
        delta = wrap_angle(phi - float(theta0))
        energy += 0.5 * float(spring) * delta * delta
    return energy


def ligand_nonbonded_energy(group, ligand_pos):
    energy = 0.0
    cutoff = float(group._v_attrs.cutoff)
    for a, b, rmin, epsilon, charge in zip(
        group.atom1[:],
        group.atom2[:],
        group.rmin[:],
        group.epsilon[:],
        group.charge[:],
    ):
        disp = ligand_pos[a] - ligand_pos[b]
        energy += charmm_pair_energy(np.dot(disp, disp), rmin, epsilon, charge, cutoff)
    return energy


def backbone_energy(group, ligand_pos, protein_pos):
    energy = 0.0
    cutoff = float(group._v_attrs.cutoff)
    protein_index = group.protein_index[:]
    ligand_rmin_half = group.ligand_rmin_half[:]
    ligand_epsilon = group.ligand_epsilon[:]
    ligand_charge = group.ligand_charge[:]
    protein_rmin_half = group.protein_rmin_half[:]
    protein_epsilon = group.protein_epsilon[:]
    protein_charge = group.protein_charge[:]
    for nl, r_lig in enumerate(ligand_pos):
        for np_idx, p_index in enumerate(protein_index):
            disp = r_lig - protein_pos[p_index, :3]
            energy += charmm_pair_energy(
                np.dot(disp, disp),
                ligand_rmin_half[nl] + protein_rmin_half[np_idx],
                math.sqrt(max(float(ligand_epsilon[nl]) * float(protein_epsilon[np_idx]), 0.0)),
                float(ligand_charge[nl]) * float(protein_charge[np_idx]),
                cutoff,
            )
    return energy


def rotamer_energy(group, ligand_pos, affine_output):
    cutoff = float(group._v_attrs.cutoff)
    affine_residue = group.affine_residue[:]
    atom_start = group.atom_start[:]
    atom_count = group.atom_count[:]
    atom_local_pos = group.atom_local_pos[:]
    atom_rmin_half = group.atom_rmin_half[:]
    atom_epsilon = group.atom_epsilon[:]
    atom_charge = group.atom_charge[:]
    ligand_rmin_half = group.ligand_rmin_half[:]
    ligand_epsilon = group.ligand_epsilon[:]
    ligand_charge = group.ligand_charge[:]

    energy = np.zeros(len(affine_residue), dtype="f8")
    for idx, residue in enumerate(affine_residue):
        count = int(atom_count[idx])
        if not count:
            continue
        start = int(atom_start[idx])
        atoms = apply_affine(affine_output[residue], atom_local_pos[start:start + count])
        e = 0.0
        for na, r_sc in enumerate(atoms, start=start):
            for nl, r_lig in enumerate(ligand_pos):
                disp = r_sc - r_lig
                e += charmm_pair_energy(
                    np.dot(disp, disp),
                    atom_rmin_half[na] + ligand_rmin_half[nl],
                    math.sqrt(max(float(atom_epsilon[na]) * float(ligand_epsilon[nl]), 0.0)),
                    float(atom_charge[na]) * float(ligand_charge[nl]),
                    cutoff,
                )
        energy[idx] = e
    return energy.reshape((-1, 1))


def snapshot_positions(handle):
    snapshots = [("initial", handle.root.input.pos[:, :, 0].astype("f4"))]
    latest = latest_output_group(handle)
    if latest is not None:
        snapshots.append(("final", latest.pos[-1, 0].astype("f4")))
    return snapshots


def compare_snapshot(handle, engine, label, pos):
    ligand_start = int(handle.root.input.ligand._v_attrs.atom_start)
    ligand_pos = pos[ligand_start:]
    protein_pos = pos[:ligand_start]
    potential = handle.root.input.potential

    engine.energy(pos.astype("f4"))
    affine_output = engine.get_output("affine_alignment")
    infer_output = engine.get_output("infer_H_O")

    print("\nSnapshot:", label)
    compare_scalar(
        "ligand_harmonic_bond",
        bond_energy(potential.ligand_harmonic_bond, ligand_pos),
        scalar_output(engine, "ligand_harmonic_bond"),
    )
    compare_scalar(
        "ligand_angle",
        angle_energy(potential.ligand_angle, ligand_pos),
        scalar_output(engine, "ligand_angle"),
    )
    compare_scalar(
        "ligand_torsion",
        torsion_energy(potential.ligand_torsion, ligand_pos),
        scalar_output(engine, "ligand_torsion"),
    )
    compare_scalar(
        "ligand_improper",
        improper_energy(potential.ligand_improper, ligand_pos),
        scalar_output(engine, "ligand_improper"),
    )
    compare_scalar(
        "charmm_ligand_nonbonded",
        ligand_nonbonded_energy(potential.charmm_ligand_nonbonded, ligand_pos),
        scalar_output(engine, "charmm_ligand_nonbonded"),
    )
    compare_scalar(
        "ligand_backbone_nonbonded",
        backbone_energy(potential.ligand_backbone_nonbonded, ligand_pos, protein_pos),
        scalar_output(engine, "ligand_backbone_nonbonded"),
    )
    compare_scalar(
        "ligand_virtual_backbone_nonbonded",
        backbone_energy(potential.ligand_virtual_backbone_nonbonded, ligand_pos, infer_output[:, :3]),
        scalar_output(engine, "ligand_virtual_backbone_nonbonded"),
    )
    compare_array(
        "ligand_rotamer_charmm_1body",
        rotamer_energy(potential.ligand_rotamer_charmm_1body, ligand_pos, affine_output),
        engine.get_output("ligand_rotamer_charmm_1body"),
    )


def main():
    print("OpenMM status:", try_openmm_status())
    ff_inputs = ci.resolve_forcefield_inputs(example_root())
    print("CHARMM topology files:")
    for path in ff_inputs["topology_paths"]:
        print("  ", path)
    print("CHARMM parameter files:")
    for path in ff_inputs["parameter_paths"]:
        print("  ", path)
    if ff_inputs["stream_paths"]:
        print("CHARMM stream files:")
        for path in ff_inputs["stream_paths"]:
            print("  ", path)
    ue = load_upside_engine()
    with tb.open_file(run_path()) as handle:
        engine = ue.Upside(run_path())
        for label, pos in snapshot_positions(handle):
            compare_snapshot(handle, engine, label, pos)
    print("\nFrozen force-field comparison passed for", run_path())


if __name__ == "__main__":
    main()
