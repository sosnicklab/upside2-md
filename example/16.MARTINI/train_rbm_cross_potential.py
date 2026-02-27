#!/usr/bin/env python3

import argparse
import json
import sys
from pathlib import Path

try:
    import h5py
    import numpy as np
except ModuleNotFoundError:
    # Allow `--help` to work even when runtime deps are not yet installed.
    if any(arg in ("-h", "--help") for arg in sys.argv[1:]):
        h5py = None
        np = None
    else:
        raise


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Train MARTINI RBM cross-potential weights from prepared .up files "
            "and optionally inject the trained node into each input file."
        )
    )
    parser.add_argument("--manifest", default=None, help="JSON manifest with training systems.")
    parser.add_argument("--up-files", nargs="*", default=[], help="Direct list of .up files.")
    parser.add_argument("--artifact-out", required=True, help="Output RBM artifact path (HDF5).")
    parser.add_argument("--epochs", type=int, default=25)
    parser.add_argument("--learning-rate", type=float, default=5e-4)
    parser.add_argument("--weight-decay", type=float, default=1e-4)
    parser.add_argument("--cutoff", type=float, default=12.0)
    parser.add_argument("--n-radial", type=int, default=12)
    parser.add_argument("--radial-min", type=float, default=2.0)
    parser.add_argument("--radial-max", type=float, default=12.0)
    parser.add_argument(
        "--radial-width",
        type=float,
        default=0.0,
        help="If <=0, use automatic spacing-based width.",
    )
    parser.add_argument("--negatives", type=int, default=1)
    parser.add_argument("--chunk-size", type=int, default=128)
    parser.add_argument("--max-protein-atoms", type=int, default=0)
    parser.add_argument("--max-env-atoms", type=int, default=0)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--inject-updates", type=int, choices=[0, 1], default=0)
    return parser.parse_args()


def decode_bytes_array(arr):
    out = []
    for item in arr:
        if isinstance(item, (bytes, np.bytes_)):
            out.append(item.decode("utf-8", errors="ignore"))
        else:
            out.append(str(item))
    return np.array(out, dtype=object)


def normalize_role(role):
    return str(role).strip().upper()


def normalize_particle_class(label):
    return str(label).strip().upper()


def protein_class_from_role(role):
    role_u = normalize_role(role)
    if role_u == "BB":
        return "BB"
    if role_u.startswith("SC"):
        return "SC"
    return "PROTEIN_OTHER"


def env_class_from_role_and_particle(role, particle_class):
    role_u = normalize_role(role)
    part_u = normalize_particle_class(particle_class)
    if part_u == "LIPID":
        if role_u:
            return f"LIPID_{role_u}"
        return "LIPID_UNK"
    if part_u:
        return part_u
    if role_u:
        return f"OTHER_{role_u}"
    return "OTHER"


def load_manifest_paths(manifest_path):
    with manifest_path.open("r", encoding="utf-8") as fh:
        data = json.load(fh)

    paths = []
    if isinstance(data, dict):
        if isinstance(data.get("systems"), list):
            items = data["systems"]
        elif isinstance(data.get("up_files"), list):
            items = data["up_files"]
        else:
            raise ValueError(f"Unsupported manifest schema: {manifest_path}")
    elif isinstance(data, list):
        items = data
    else:
        raise ValueError(f"Unsupported manifest top-level type: {type(data)}")

    for item in items:
        if isinstance(item, str):
            paths.append(Path(item).expanduser().resolve())
            continue
        if isinstance(item, dict):
            up_file = item.get("up_file") or item.get("upside_input")
            if not up_file:
                continue
            paths.append(Path(up_file).expanduser().resolve())
            continue
    return paths


def collect_up_paths(args):
    paths = []

    if args.manifest:
        manifest_path = Path(args.manifest).expanduser().resolve()
        if not manifest_path.exists():
            raise FileNotFoundError(f"Manifest not found: {manifest_path}")
        paths.extend(load_manifest_paths(manifest_path))

    for raw in args.up_files:
        paths.append(Path(raw).expanduser().resolve())

    # Stable unique ordering.
    uniq = []
    seen = set()
    for p in paths:
        if p in seen:
            continue
        seen.add(p)
        uniq.append(p)

    if not uniq:
        raise ValueError("No input .up files provided. Use --manifest and/or --up-files.")

    missing = [p for p in uniq if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing .up files:\n" + "\n".join(str(p) for p in missing))

    return uniq


def read_positions(h5):
    pos = h5["/input/pos"][:]
    if pos.ndim == 3:
        # Expected shape (n_atom, 3, 1)
        return pos[:, :, 0].astype(np.float64, copy=False)
    if pos.ndim == 2:
        return pos.astype(np.float64, copy=False)
    raise ValueError(f"Unsupported /input/pos shape: {pos.shape}")


def extract_system_data(up_path, rng, max_protein_atoms=0, max_env_atoms=0):
    with h5py.File(up_path, "r") as h5:
        xyz = read_positions(h5)
        n_atom = xyz.shape[0]

        if "/input/atom_roles" in h5:
            atom_roles = decode_bytes_array(h5["/input/atom_roles"][:])
        elif "/input/atom_names" in h5:
            atom_roles = decode_bytes_array(h5["/input/atom_names"][:])
        else:
            atom_roles = np.array(["UNK"] * n_atom, dtype=object)

        if "/input/particle_class" in h5:
            particle_class = decode_bytes_array(h5["/input/particle_class"][:])
        else:
            particle_class = np.array(["OTHER"] * n_atom, dtype=object)

        if "/input/hybrid_env_topology/protein_membership" in h5:
            membership = h5["/input/hybrid_env_topology/protein_membership"][:].astype(np.int32)
            if membership.shape[0] != n_atom:
                raise ValueError(
                    f"{up_path}: protein_membership length {membership.shape[0]} != n_atom {n_atom}"
                )
            protein_mask = membership >= 0
        else:
            protein_mask = np.array(
                [normalize_particle_class(c) == "PROTEIN" for c in particle_class], dtype=bool
            )

        env_mask = ~protein_mask
        protein_indices = np.where(protein_mask)[0].astype(np.int32)
        env_indices = np.where(env_mask)[0].astype(np.int32)

        if max_protein_atoms > 0 and protein_indices.size > max_protein_atoms:
            sel = rng.choice(protein_indices.size, size=max_protein_atoms, replace=False)
            protein_indices = np.sort(protein_indices[sel])
        if max_env_atoms > 0 and env_indices.size > max_env_atoms:
            sel = rng.choice(env_indices.size, size=max_env_atoms, replace=False)
            env_indices = np.sort(env_indices[sel])

        if protein_indices.size == 0:
            raise ValueError(f"{up_path}: no protein atoms identified")
        if env_indices.size == 0:
            raise ValueError(f"{up_path}: no environment atoms identified")

        protein_labels = [protein_class_from_role(atom_roles[i]) for i in protein_indices]
        env_labels = [
            env_class_from_role_and_particle(atom_roles[i], particle_class[i]) for i in env_indices
        ]

        return {
            "up_path": up_path,
            "xyz": xyz,
            "protein_atom_indices": protein_indices,
            "env_atom_indices": env_indices,
            "protein_labels": protein_labels,
            "env_labels": env_labels,
        }


def build_radial_basis(args):
    if args.n_radial <= 0:
        raise ValueError("--n-radial must be > 0")
    if args.radial_max <= args.radial_min:
        raise ValueError("--radial-max must be > --radial-min")

    centers = np.linspace(args.radial_min, args.radial_max, args.n_radial, dtype=np.float64)
    if args.radial_width > 0:
        width = float(args.radial_width)
    elif args.n_radial > 1:
        width = float((args.radial_max - args.radial_min) / (args.n_radial - 1))
    else:
        width = float(max(1.0, args.radial_max - args.radial_min))

    widths = np.full(args.n_radial, width, dtype=np.float64)
    return centers, widths


def accumulate_features(
    protein_xyz,
    protein_class_idx,
    env_xyz,
    env_class_idx,
    n_protein_class,
    n_env_class,
    centers,
    widths,
    cutoff,
    chunk_size,
):
    n_radial = centers.shape[0]
    feat = np.zeros((n_protein_class, n_env_class, n_radial), dtype=np.float64)
    cutoff_sq = float(cutoff) * float(cutoff)

    protein_by_class = [np.where(protein_class_idx == c)[0] for c in range(n_protein_class)]
    env_by_class = [np.where(env_class_idx == c)[0] for c in range(n_env_class)]

    for pc, psel in enumerate(protein_by_class):
        if psel.size == 0:
            continue
        for ec, esel in enumerate(env_by_class):
            if esel.size == 0:
                continue

            env_block = env_xyz[esel]
            for start in range(0, psel.size, chunk_size):
                pblock = protein_xyz[psel[start : start + chunk_size]]
                diff = pblock[:, None, :] - env_block[None, :, :]
                r2 = np.einsum("ijk,ijk->ij", diff, diff)
                mask = (r2 < cutoff_sq) & (r2 > 1e-8)
                if not np.any(mask):
                    continue

                r = np.sqrt(np.maximum(r2, 0.0))
                delta = (r[:, :, None] - centers[None, None, :]) / widths[None, None, :]
                basis = np.exp(-0.5 * delta * delta)
                basis *= mask[:, :, None]
                feat[pc, ec, :] += np.sum(basis, axis=(0, 1))

    return feat


def write_artifact(artifact_path, args, centers, widths, protein_classes, env_classes, weights):
    artifact_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(artifact_path, "w") as h5:
        meta = h5.create_group("meta")
        meta.attrs["schema"] = "martini_rbm_cross_v1"
        meta.attrs["epochs"] = np.int32(args.epochs)
        meta.attrs["learning_rate"] = np.float32(args.learning_rate)
        meta.attrs["weight_decay"] = np.float32(args.weight_decay)
        meta.attrs["cutoff"] = np.float32(args.cutoff)
        meta.attrs["negatives"] = np.int32(args.negatives)

        classes = h5.create_group("classes")
        sdt = h5py.string_dtype(encoding="utf-8")
        classes.create_dataset("protein", data=np.array(protein_classes, dtype=sdt))
        classes.create_dataset("environment", data=np.array(env_classes, dtype=sdt))

        rbm = h5.create_group("rbm")
        cross = rbm.create_group("cross")
        cross.attrs["cutoff"] = np.float32(args.cutoff)
        cross.create_dataset("radial_centers", data=centers.astype(np.float32))
        cross.create_dataset("radial_widths", data=widths.astype(np.float32))
        cross.create_dataset("weights", data=weights.astype(np.float32))

        # Compatibility views expected by current hybrid tooling.
        for name in ("bb", "sc"):
            grp = rbm.create_group(name)
            grp.create_dataset("radial_centers", data=centers.astype(np.float32))
            grp.create_dataset("radial_widths", data=widths.astype(np.float32))
            grp.create_dataset("weights", data=weights.astype(np.float32))

        views = h5.create_group("views")
        views.create_dataset("protein_classes", data=np.array(protein_classes, dtype=sdt))
        views.create_dataset("environment_classes", data=np.array(env_classes, dtype=sdt))

        priors = h5.create_group("priors")
        priors.create_dataset("weight_l2", data=np.array([args.weight_decay], dtype=np.float32))


def inject_node_into_up(
    up_path,
    cutoff,
    centers,
    widths,
    weights,
    protein_atom_indices,
    protein_class_indices,
    env_atom_indices,
    env_class_indices,
):
    with h5py.File(up_path, "r+") as h5:
        input_grp = h5.require_group("input")
        pot_grp = input_grp.require_group("potential")
        node_name = "martini_rbm_cross_potential"
        if node_name in pot_grp:
            del pot_grp[node_name]
        grp = pot_grp.create_group(node_name)

        grp.attrs["arguments"] = np.array([b"pos"])
        grp.attrs["initialized"] = np.int8(1)
        grp.attrs["cutoff"] = np.float32(cutoff)

        grp.create_dataset("protein_atom_indices", data=protein_atom_indices.astype(np.int32))
        grp.create_dataset("protein_class_index", data=protein_class_indices.astype(np.int32))
        grp.create_dataset("env_atom_indices", data=env_atom_indices.astype(np.int32))
        grp.create_dataset("env_class_index", data=env_class_indices.astype(np.int32))
        grp.create_dataset("radial_centers", data=centers.astype(np.float32))
        grp.create_dataset("radial_widths", data=widths.astype(np.float32))
        grp.create_dataset("weights", data=weights.astype(np.float32))


def main():
    args = parse_args()
    if h5py is None or np is None:
        raise ModuleNotFoundError(
            "train_rbm_cross_potential.py requires numpy and h5py at runtime. "
            "Install dependencies in the active environment."
        )
    if args.cutoff <= 0:
        raise ValueError("--cutoff must be > 0")
    if args.negatives <= 0:
        raise ValueError("--negatives must be > 0")
    if args.chunk_size <= 0:
        raise ValueError("--chunk-size must be > 0")

    rng = np.random.default_rng(args.seed)
    up_paths = collect_up_paths(args)

    systems = []
    for up_path in up_paths:
        systems.append(
            extract_system_data(
                up_path,
                rng=rng,
                max_protein_atoms=max(0, args.max_protein_atoms),
                max_env_atoms=max(0, args.max_env_atoms),
            )
        )

    protein_classes = sorted({label for s in systems for label in s["protein_labels"]})
    env_classes = sorted({label for s in systems for label in s["env_labels"]})
    protein_class_to_idx = {name: i for i, name in enumerate(protein_classes)}
    env_class_to_idx = {name: i for i, name in enumerate(env_classes)}

    n_protein_class = len(protein_classes)
    n_env_class = len(env_classes)

    for s in systems:
        s["protein_class_index"] = np.array(
            [protein_class_to_idx[name] for name in s["protein_labels"]], dtype=np.int32
        )
        s["env_class_index"] = np.array(
            [env_class_to_idx[name] for name in s["env_labels"]], dtype=np.int32
        )
        s["protein_xyz"] = s["xyz"][s["protein_atom_indices"]]
        s["env_xyz"] = s["xyz"][s["env_atom_indices"]]

    centers, widths = build_radial_basis(args)
    n_radial = centers.shape[0]
    weights = np.zeros((n_protein_class, n_env_class, n_radial), dtype=np.float64)

    for epoch in range(args.epochs):
        diff_accum = np.zeros_like(weights)
        delta_accum = 0.0

        for s in systems:
            f_pos = accumulate_features(
                s["protein_xyz"],
                s["protein_class_index"],
                s["env_xyz"],
                s["env_class_index"],
                n_protein_class,
                n_env_class,
                centers,
                widths,
                args.cutoff,
                args.chunk_size,
            )

            f_neg = np.zeros_like(f_pos)
            for _ in range(args.negatives):
                shuffled_env = s["env_xyz"][rng.permutation(s["env_xyz"].shape[0])]
                f_neg += accumulate_features(
                    s["protein_xyz"],
                    s["protein_class_index"],
                    shuffled_env,
                    s["env_class_index"],
                    n_protein_class,
                    n_env_class,
                    centers,
                    widths,
                    args.cutoff,
                    args.chunk_size,
                )
            f_neg /= float(args.negatives)

            diff = f_pos - f_neg
            diff_accum += diff
            delta_accum += float(np.sum(weights * diff))

        batch_scale = 1.0 / float(len(systems))
        grad = args.weight_decay * weights - diff_accum * batch_scale
        weights -= args.learning_rate * grad

        loss = 0.5 * args.weight_decay * float(np.sum(weights * weights)) - delta_accum * batch_scale
        print(
            f"epoch={epoch + 1:03d}/{args.epochs:03d} "
            f"loss={loss:.6f} avg_delta={delta_accum * batch_scale:.6f}"
        )

    artifact_path = Path(args.artifact_out).expanduser().resolve()
    write_artifact(
        artifact_path,
        args=args,
        centers=centers,
        widths=widths,
        protein_classes=protein_classes,
        env_classes=env_classes,
        weights=weights,
    )
    print(f"Wrote RBM artifact: {artifact_path}")

    if args.inject_updates:
        for s in systems:
            inject_node_into_up(
                up_path=s["up_path"],
                cutoff=args.cutoff,
                centers=centers,
                widths=widths,
                weights=weights,
                protein_atom_indices=s["protein_atom_indices"],
                protein_class_indices=s["protein_class_index"],
                env_atom_indices=s["env_atom_indices"],
                env_class_indices=s["env_class_index"],
            )
            print(f"Injected martini_rbm_cross_potential into: {s['up_path']}")


if __name__ == "__main__":
    main()
