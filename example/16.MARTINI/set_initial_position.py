#!/usr/bin/env python
import sys
import os
import h5py
import numpy as np

TWOPI = 2.0 * np.pi


def decode_stage_label(value):
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="ignore").strip().lower()
    return str(value).strip().lower()


def env_bool(name, default=False):
    raw = os.getenv(name)
    if raw is None:
        return default
    value = str(raw).strip().lower()
    return value not in ("0", "false", "no", "off", "")


def recenter_protein_for_production(h5f, pos, box_lengths):
    """Recenters protein BB to box center in production handoff."""
    if box_lengths is None:
        return pos
    if '/input/stage_parameters' not in h5f:
        return pos
    stage = decode_stage_label(h5f['/input/stage_parameters'].attrs.get('current_stage', b''))
    if stage != 'production':
        return pos
    if '/input/hybrid_bb_map/bb_atom_index' not in h5f:
        return pos

    L = np.asarray(box_lengths, dtype=np.float64).reshape(3)
    if not np.all(np.isfinite(L)) or np.any(L <= 0.0):
        return pos

    bb_idx = h5f['/input/hybrid_bb_map/bb_atom_index'][:].astype(np.int32)
    bb_idx = bb_idx[(bb_idx >= 0) & (bb_idx < pos.shape[0])]
    if bb_idx.size == 0:
        return pos

    bb_xyz = pos[bb_idx, :, 0].astype(np.float64)
    bb_center = np.zeros(3, dtype=np.float64)
    for ax in range(3):
        coords = np.mod(bb_xyz[:, ax], L[ax])
        angles = TWOPI * coords / L[ax]
        s = np.mean(np.sin(angles))
        c = np.mean(np.cos(angles))
        if abs(s) < 1e-12 and abs(c) < 1e-12:
            bb_center[ax] = float(np.mean(coords))
        else:
            ang = np.arctan2(s, c)
            if ang < 0.0:
                ang += TWOPI
            bb_center[ax] = L[ax] * ang / TWOPI

    target_center = 0.5 * L
    delta = target_center - bb_center
    # Apply a rigid translation to all atoms without per-atom modulo wrapping.
    # Per-atom wrapping can alter raw pair displacements in code paths that do
    # not apply minimum-image, causing artificial energy jumps at stage handoff.
    pos[:, :, 0] = pos[:, :, 0] + delta.astype(pos.dtype, copy=False)[None, :]

    print(
        "Recentered production protein BB to box center: "
        f"before=({bb_center[0]:.3f}, {bb_center[1]:.3f}, {bb_center[2]:.3f}) "
        f"target=({target_center[0]:.3f}, {target_center[1]:.3f}, {target_center[2]:.3f})"
    )
    return pos


def refresh_hybrid_reference_carriers(h5f, pos):
    if '/input/hybrid_bb_map' not in h5f:
        return pos

    grp = h5f['/input/hybrid_bb_map']
    required = ('bb_atom_index', 'atom_indices', 'weights', 'reference_atom_coords')
    if not all((f'/input/hybrid_bb_map/{k}' in h5f) for k in required):
        return pos

    bb_idx = grp['bb_atom_index'][:].astype(np.int32)
    comp_idx = grp['atom_indices'][:].astype(np.int32)
    comp_w = grp['weights'][:].astype(np.float64)
    ref_xyz = grp['reference_atom_coords'][:].astype(np.float64)

    if comp_idx.ndim != 2 or comp_idx.shape[1] != 4:
        return pos
    if ref_xyz.shape != (comp_idx.shape[0], 4, 3):
        return pos

    n_atom = pos.shape[0]
    for k in range(comp_idx.shape[0]):
        bb = int(bb_idx[k]) if k < bb_idx.shape[0] else -1
        if bb < 0 or bb >= n_atom:
            continue

        valid = (comp_idx[k] >= 0) & (comp_idx[k] < n_atom) & (comp_w[k] > 0.0)
        if not np.any(valid):
            continue

        w = comp_w[k][valid]
        wsum = float(np.sum(w))
        if wsum <= 0.0:
            continue
        w = w / wsum

        ref_pts = ref_xyz[k][valid]
        ref_com = np.sum(ref_pts * w[:, None], axis=0)
        bb_pos = pos[bb, :, 0].astype(np.float64)
        shift = bb_pos - ref_com

        target_idx = comp_idx[k][valid]
        aligned = ref_pts + shift[None, :]
        for ai, pxyz in zip(target_idx, aligned):
            pos[int(ai), :, 0] = pxyz.astype(pos.dtype)

    return pos

def set_initial_position(input_file, output_file):
    strict_copy = env_bool("UPSIDE_SET_INITIAL_STRICT_COPY", False)
    apply_refresh_hybrid = env_bool(
        "UPSIDE_SET_INITIAL_REFRESH_HYBRID_CARRIERS",
        default=(not strict_copy),
    )
    apply_recenter_production = env_bool(
        "UPSIDE_SET_INITIAL_RECENTER_PRODUCTION",
        default=(not strict_copy),
    )

    # Open input file and get last frame's position and box dimensions
    with h5py.File(input_file, 'r') as f:
        # Try to read from output/pos (if simulation has already produced frames)
        if '/output/pos' in f and f['/output/pos'].shape[0] > 0:
            last_pos = f['/output/pos'][-1, 0, :, :]  # Last frame, first replica
            last_pos = last_pos[:, :, np.newaxis]  # Add replica dimension (1)
        else:
            # Fallback when output exists but is empty (e.g. minimize-only runs)
            last_pos = f['/input/pos'][:, :, 0]
            last_pos = last_pos[:, :, np.newaxis]

        # Get last frame's box dimensions if available
        last_box = None
        if '/output/box' in f:
            box_data = f['/output/box'][:]
            if box_data.size > 0:
                last_box = box_data[-1]
                # Handle different box data formats
                if len(last_box.shape) == 2 and last_box.shape[1] == 3:
                    last_box = last_box[0]
        # If no output/box, try to get box from input potential attributes
        if last_box is None:
            if '/input/potential/martini_potential' in f:
                pot_grp = f['/input/potential/martini_potential']
                if all(k in pot_grp.attrs for k in ('x_len', 'y_len', 'z_len')):
                    last_box = np.array([
                        pot_grp.attrs['x_len'],
                        pot_grp.attrs['y_len'],
                        pot_grp.attrs['z_len']
                    ])

    # Open output file and write last frame's position to input/pos
    with h5py.File(output_file, 'r+') as f:
        target_pos = f['/input/pos'][:]
        target_n = target_pos.shape[0]
        source_n = last_pos.shape[0]
        if source_n != target_n:
            merged = target_pos.copy()
            n_copy = min(source_n, target_n)
            merged[:n_copy, :, :] = last_pos[:n_copy, :, :]
            last_pos = merged

        if strict_copy and not apply_refresh_hybrid and not apply_recenter_production:
            print("Strict handoff mode: preserving exact coordinates from previous stage output.")

        # If hybrid reference carriers are present in the target file, align them
        # to the current BB proxy positions after stage handoff.
        if apply_refresh_hybrid:
            last_pos = refresh_hybrid_reference_carriers(f, last_pos)
        if apply_recenter_production:
            last_pos = recenter_protein_for_production(f, last_pos, last_box)

        if '/input/pos' in f:
            del f['/input/pos']
        f.create_dataset('/input/pos', data=last_pos)

        # Update box dimensions in martini_potential if available
        if last_box is not None and '/input/potential/martini_potential' in f:
            pot_grp = f['/input/potential/martini_potential']
            pot_grp.attrs['x_len'] = float(last_box[0])
            pot_grp.attrs['y_len'] = float(last_box[1])
            pot_grp.attrs['z_len'] = float(last_box[2])
            print(f"Updated box dimensions: x={last_box[0]:.3f}, y={last_box[1]:.3f}, z={last_box[2]:.3f}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 set_initial_position.py <input_file> <output_file>")
        sys.exit(1)
    set_initial_position(sys.argv[1], sys.argv[2])
