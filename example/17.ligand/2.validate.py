#!/usr/bin/env python3

import os
import sys

import numpy as np
import tables as tb


SIM_ID = "smoke"
MAX_ANCHOR_RMS = 0.5
MAX_POCKET_ALIGNED_LIGAND_RMSD = 3.0


def kabsch_row(mobile, target):
    center_mobile = mobile.mean(axis=0)
    center_target = target.mean(axis=0)
    mobile_centered = mobile - center_mobile
    target_centered = target - center_target
    cov = mobile_centered.T.dot(target_centered)
    u, _s, vt = np.linalg.svd(cov)
    rot = u.dot(vt)
    if np.linalg.det(rot) < 0.0:
        u[:, -1] *= -1.0
        rot = u.dot(vt)
    return rot, center_mobile, center_target


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
    if latest is None:
        raise RuntimeError("No output group found in %s" % run_path())
    return latest


def main():
    ue = load_upside_engine()
    with tb.open_file(run_path()) as t:
        ligand_start = int(t.root.input.ligand._v_attrs.atom_start)
        ligand_initial = t.root.input.ligand.initial_pos[:]
        protein_initial = t.root.input.pos[:ligand_start, :, 0]
        pocket_residue = t.root.input.ligand.pocket_residue[:]
        potential_names = set(t.root.input.potential._v_children)
        anchor_args = [x.decode("ascii") for x in t.root.input.potential.Distance3D_ligand_anchor._v_attrs.arguments]
        anchor_id = t.root.input.potential.Distance3D_ligand_anchor.id[:]
        anchor_eq = t.root.input.potential.Spring_ligand_anchor.equil_dist[:]
        out = latest_output_group(t)
        pos = out.pos[:, 0]

        if pos.shape[0] < 2:
            raise RuntimeError("Expected multiple frames, found %d" % pos.shape[0])
        if not np.isfinite(pos).all():
            raise RuntimeError("Trajectory contains non-finite coordinates")

        ligand_xyz = pos[:, ligand_start:]
        if ligand_xyz.shape[1] != ligand_initial.shape[0]:
            raise RuntimeError("Ligand atom count mismatch in trajectory")
        required_nodes = {
            "charmm_ligand_nonbonded",
            "ligand_backbone_nonbonded",
            "ligand_virtual_backbone_nonbonded",
            "ligand_improper",
            "ligand_rotamer_charmm_1body",
            "infer_H_O",
        }
        if not required_nodes.issubset(potential_names):
            raise RuntimeError("Physical ligand nodes missing from prepared config")
        rotamer_datasets = set(t.root.input.potential.ligand_rotamer_charmm_1body._v_children)
        required_rotamer_datasets = {
            "atom_start",
            "atom_count",
            "atom_name",
            "atom_element",
            "atom_local_pos",
            "atom_rmin_half",
            "atom_epsilon",
            "atom_charge",
        }
        if not required_rotamer_datasets.issubset(rotamer_datasets):
            raise RuntimeError("Atomistic rotamer sidechain library missing from prepared config")
        if anchor_args != ["slice_ligand", "placement_fixed_point_vector_only"]:
            raise RuntimeError("Ligand anchors are not sidechain-placement based")

        protein_xyz = pos[:, :ligand_start]
        pocket_atoms = np.concatenate(
            [np.arange(residue * 3, residue * 3 + 3, dtype="i4") for residue in pocket_residue]
        )
        rot, center_mobile, center_target = kabsch_row(protein_xyz[-1, pocket_atoms], protein_initial[pocket_atoms])
        ligand_final_aligned = (ligand_xyz[-1] - center_mobile[None, :]).dot(rot) + center_target[None, :]
        protein_final_aligned = (protein_xyz[-1] - center_mobile[None, :]).dot(rot) + center_target[None, :]
        ligand_rmsd = np.sqrt(np.mean(np.sum((ligand_final_aligned - ligand_initial) ** 2, axis=1)))
        pocket_rmsd = np.sqrt(
            np.mean(np.sum((protein_final_aligned[pocket_atoms] - protein_initial[pocket_atoms]) ** 2, axis=1))
        )

        engine = ue.Upside(run_path())
        engine.energy(pos[-1].astype("f4"))
        placement_pos = engine.get_output("placement_fixed_point_vector_only")[:, :3]
        deviations = []
        for ligand_atom, placement_atom in anchor_id:
            disp = pos[-1, ligand_start + ligand_atom] - placement_pos[placement_atom]
            deviations.append(np.sqrt((disp * disp).sum()))
        deviations = np.asarray(deviations, dtype="f4") - anchor_eq
        anchor_rms = np.sqrt(np.mean(deviations * deviations))
        print("Ligand anchor source:", anchor_args[1])
        print("Final pocket-aligned ligand RMSD (A): %.3f" % ligand_rmsd)
        print("Final pocket backbone RMSD (A): %.3f" % pocket_rmsd)
        print("Final ligand anchor distance RMS error (A): %.3f" % anchor_rms)
        if ligand_rmsd > MAX_POCKET_ALIGNED_LIGAND_RMSD:
            raise RuntimeError("Ligand pose drifted too far after pocket alignment")
        if anchor_rms > MAX_ANCHOR_RMS:
            raise RuntimeError("Ligand anchor geometry drifted too far from the starting bound pose")

    print("Validation passed for", run_path())


if __name__ == "__main__":
    main()
