import math

import numpy as np


DEG = math.pi / 180.0

RESTYPE_ORDER = (
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
)


def canonical_restype(restype):
    return "PRO" if restype == "CPR" else restype


def _param(name, rmin_half, epsilon, charge):
    return {
        "name": name,
        "element": name[0],
        "rmin_half": float(rmin_half),
        "epsilon": float(epsilon),
        "charge": float(charge),
    }


RESIDUE_ATOMS = {
    "ALA": [_param("CB", 2.00, 0.07, 0.00)],
    "ARG": [
        _param("CB", 2.00, 0.07, 0.00), _param("CG", 2.00, 0.07, 0.00),
        _param("CD", 2.00, 0.07, 0.00), _param("NE", 1.85, 0.20, 0.00),
        _param("CZ", 2.00, 0.07, 0.50), _param("NH1", 1.85, 0.20, 0.25),
        _param("NH2", 1.85, 0.20, 0.25),
    ],
    "ASN": [
        _param("CB", 2.00, 0.07, 0.00), _param("CG", 2.00, 0.07, 0.40),
        _param("OD1", 1.70, 0.12, -0.40), _param("ND2", 1.85, 0.20, 0.00),
    ],
    "ASP": [
        _param("CB", 2.00, 0.07, 0.00), _param("CG", 2.00, 0.07, 0.00),
        _param("OD1", 1.70, 0.12, -0.50), _param("OD2", 1.70, 0.12, -0.50),
    ],
    "CYS": [_param("CB", 2.00, 0.07, 0.00), _param("SG", 2.00, 0.25, 0.00)],
    "GLN": [
        _param("CB", 2.00, 0.07, 0.00), _param("CG", 2.00, 0.07, 0.00),
        _param("CD", 2.00, 0.07, 0.40), _param("OE1", 1.70, 0.12, -0.40),
        _param("NE2", 1.85, 0.20, 0.00),
    ],
    "GLU": [
        _param("CB", 2.00, 0.07, 0.00), _param("CG", 2.00, 0.07, 0.00),
        _param("CD", 2.00, 0.07, 0.00), _param("OE1", 1.70, 0.12, -0.50),
        _param("OE2", 1.70, 0.12, -0.50),
    ],
    "GLY": [],
    "HIS": [
        _param("CB", 2.00, 0.07, 0.00), _param("CG", 1.99, 0.07, 0.00),
        _param("ND1", 1.85, 0.20, 0.00), _param("CD2", 1.99, 0.07, 0.00),
        _param("CE1", 1.99, 0.07, 0.00), _param("NE2", 1.85, 0.20, 0.00),
    ],
    "ILE": [
        _param("CB", 2.00, 0.07, 0.00), _param("CG1", 2.00, 0.07, 0.00),
        _param("CG2", 2.00, 0.07, 0.00), _param("CD1", 2.00, 0.07, 0.00),
    ],
    "LEU": [
        _param("CB", 2.00, 0.07, 0.00), _param("CG", 2.00, 0.07, 0.00),
        _param("CD1", 2.00, 0.07, 0.00), _param("CD2", 2.00, 0.07, 0.00),
    ],
    "LYS": [
        _param("CB", 2.00, 0.07, 0.00), _param("CG", 2.00, 0.07, 0.00),
        _param("CD", 2.00, 0.07, 0.00), _param("CE", 2.00, 0.07, 0.00),
        _param("NZ", 1.85, 0.20, 1.00),
    ],
    "MET": [
        _param("CB", 2.00, 0.07, 0.00), _param("CG", 2.00, 0.07, 0.00),
        _param("SD", 2.00, 0.25, 0.00), _param("CE", 2.00, 0.07, 0.00),
    ],
    "PHE": [
        _param("CB", 2.00, 0.07, 0.00), _param("CG", 1.99, 0.07, 0.00),
        _param("CD1", 1.99, 0.07, 0.00), _param("CD2", 1.99, 0.07, 0.00),
        _param("CE1", 1.99, 0.07, 0.00), _param("CE2", 1.99, 0.07, 0.00),
        _param("CZ", 1.99, 0.07, 0.00),
    ],
    "PRO": [
        _param("CB", 2.00, 0.07, 0.00), _param("CG", 2.00, 0.07, 0.00),
        _param("CD", 2.00, 0.07, 0.00),
    ],
    "SER": [_param("CB", 2.00, 0.07, 0.00), _param("OG", 1.70, 0.12, 0.00)],
    "THR": [
        _param("CB", 2.00, 0.07, 0.00), _param("OG1", 1.70, 0.12, 0.00),
        _param("CG2", 2.00, 0.07, 0.00),
    ],
    "TRP": [
        _param("CB", 2.00, 0.07, 0.00), _param("CG", 1.99, 0.07, 0.00),
        _param("CD1", 1.99, 0.07, 0.00), _param("CD2", 1.99, 0.07, 0.00),
        _param("NE1", 1.85, 0.20, 0.00), _param("CE2", 1.99, 0.07, 0.00),
        _param("CE3", 1.99, 0.07, 0.00), _param("CZ2", 1.99, 0.07, 0.00),
        _param("CZ3", 1.99, 0.07, 0.00), _param("CH2", 1.99, 0.07, 0.00),
    ],
    "TYR": [
        _param("CB", 2.00, 0.07, 0.00), _param("CG", 1.99, 0.07, 0.00),
        _param("CD1", 1.99, 0.07, 0.00), _param("CD2", 1.99, 0.07, 0.00),
        _param("CE1", 1.99, 0.07, 0.00), _param("CE2", 1.99, 0.07, 0.00),
        _param("CZ", 1.99, 0.07, 0.00), _param("OH", 1.70, 0.12, 0.00),
    ],
    "VAL": [
        _param("CB", 2.00, 0.07, 0.00), _param("CG1", 2.00, 0.07, 0.00),
        _param("CG2", 2.00, 0.07, 0.00),
    ],
}

RESIDUE_BONDS = {
    "ALA": [],
    "ARG": [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (4, 6)],
    "ASN": [(0, 1), (1, 2), (1, 3)],
    "ASP": [(0, 1), (1, 2), (1, 3)],
    "CYS": [(0, 1)],
    "GLN": [(0, 1), (1, 2), (2, 3), (2, 4)],
    "GLU": [(0, 1), (1, 2), (2, 3), (2, 4)],
    "GLY": [],
    "HIS": [(0, 1), (1, 2), (1, 3), (2, 4), (3, 5), (4, 5)],
    "ILE": [(0, 1), (0, 2), (1, 3)],
    "LEU": [(0, 1), (1, 2), (1, 3)],
    "LYS": [(0, 1), (1, 2), (2, 3), (3, 4)],
    "MET": [(0, 1), (1, 2), (2, 3)],
    "PHE": [(0, 1), (1, 2), (1, 3), (2, 4), (3, 5), (4, 6), (5, 6)],
    "PRO": [(0, 1), (1, 2)],
    "SER": [(0, 1)],
    "THR": [(0, 1), (0, 2)],
    "TRP": [(0, 1), (1, 2), (1, 3), (2, 4), (3, 5), (3, 6), (4, 5), (5, 7), (6, 8), (7, 9), (8, 9)],
    "TYR": [(0, 1), (1, 2), (1, 3), (2, 4), (3, 5), (4, 6), (5, 6), (6, 7)],
    "VAL": [(0, 1), (0, 2)],
}


def fill_missing_chi(restype, chi_in):
    chi = np.asarray(chi_in, dtype="f4").copy()
    trans = math.pi
    if restype in ("GLN", "GLU"):
        chi[2] = 0.0
    elif restype == "LYS":
        chi[2] = trans
        chi[3] = trans
    elif restype == "ARG":
        chi[2] = trans
        chi[3] = 0.0
    elif restype == "MET":
        chi[2] = trans
    return chi


def _make_tab(phi, theta, bond):
    cp = math.cos(phi)
    sp = math.sin(phi)
    ct = math.cos(theta)
    st = math.sin(theta)
    out = np.eye(4, dtype="f4")
    out[0, 0] = -ct
    out[0, 1] = -st
    out[0, 3] = -bond * ct
    out[1, 0] = cp * st
    out[1, 1] = -cp * ct
    out[1, 2] = -sp
    out[1, 3] = bond * cp * st
    out[2, 0] = sp * st
    out[2, 1] = -sp * ct
    out[2, 2] = cp
    out[2, 3] = bond * sp * st
    return out


def _place_bb(psi, include_cb=True):
    a = np.eye(4, dtype="f4")
    a[:3, :4] = np.array([
        [0.8191292, -0.3103239, 0.4824173, -1.2079210],
        [0.5736088, 0.4423396, -0.6894263, -0.2636016],
        [0.0005532, 0.8414480, 0.5403378, -0.0009170],
    ], dtype="f4")

    n = a @ _make_tab(0.0, 0.0, 0.0)
    ca = n @ _make_tab(0.0, 0.0, 1.45)
    c = ca @ _make_tab(122.7 * DEG, 110.3 * DEG, 1.53)
    o = c @ _make_tab(psi + 180.0 * DEG, 120.5 * DEG, 1.23)
    cb = ca @ _make_tab(0.0, 110.6 * DEG, 1.53)

    pos = [n[:3, 3], ca[:3, 3], c[:3, 3], o[:3, 3]]
    if include_cb:
        pos.append(cb[:3, 3])
    return np.asarray(pos, dtype="f4"), cb


def _append(pos_list, transform):
    pos_list.append(transform[:3, 3].astype("f4"))


def _build_full_positions(restype, chi, psi=0.0):
    restype = canonical_restype(restype)
    chi = fill_missing_chi(restype, chi)

    if restype == "GLY":
        pos, _cb = _place_bb(psi, include_cb=False)
        return pos

    pos, cb = _place_bb(psi, include_cb=True)
    atoms = [row.copy() for row in pos]

    if restype == "ALA":
        return np.asarray(atoms, dtype="f4")
    if restype == "ARG":
        cg = cb @ _make_tab(chi[0], 113.9 * DEG, 1.52); _append(atoms, cg)
        cd = cg @ _make_tab(chi[1], 111.7 * DEG, 1.52); _append(atoms, cd)
        ne = cd @ _make_tab(chi[2], 111.7 * DEG, 1.46); _append(atoms, ne)
        cz = ne @ _make_tab(chi[3], 124.7 * DEG, 1.33); _append(atoms, cz)
        nh1 = cz @ _make_tab(0.0 * DEG, 120.7 * DEG, 1.33); _append(atoms, nh1)
        nh2 = cz @ _make_tab(-180.0 * DEG, 119.6 * DEG, 1.33); _append(atoms, nh2)
    elif restype == "ASN":
        cg = cb @ _make_tab(chi[0], 112.7 * DEG, 1.52); _append(atoms, cg)
        od1 = cg @ _make_tab(chi[1], 120.9 * DEG, 1.23); _append(atoms, od1)
        nd2 = cg @ _make_tab(chi[1] + 180.0 * DEG, 116.5 * DEG, 1.33); _append(atoms, nd2)
    elif restype == "ASP":
        cg = cb @ _make_tab(chi[0], 113.0 * DEG, 1.52); _append(atoms, cg)
        od1 = cg @ _make_tab(chi[1], 119.2 * DEG, 1.25); _append(atoms, od1)
        od2 = cg @ _make_tab(chi[1] - 179.9 * DEG, 118.2 * DEG, 1.25); _append(atoms, od2)
    elif restype == "CYS":
        sg = cb @ _make_tab(chi[0], 113.8 * DEG, 1.81); _append(atoms, sg)
    elif restype == "GLN":
        cg = cb @ _make_tab(chi[0], 113.9 * DEG, 1.52); _append(atoms, cg)
        cd = cg @ _make_tab(chi[1], 112.8 * DEG, 1.52); _append(atoms, cd)
        oe1 = cd @ _make_tab(chi[2], 120.9 * DEG, 1.23); _append(atoms, oe1)
        ne2 = cd @ _make_tab(chi[2] - 180.0 * DEG, 116.5 * DEG, 1.33); _append(atoms, ne2)
    elif restype == "GLU":
        cg = cb @ _make_tab(chi[0], 113.9 * DEG, 1.52); _append(atoms, cg)
        cd = cg @ _make_tab(chi[1], 113.2 * DEG, 1.52); _append(atoms, cd)
        oe1 = cd @ _make_tab(chi[2], 119.0 * DEG, 1.25); _append(atoms, oe1)
        oe2 = cd @ _make_tab(chi[2] - 180.0 * DEG, 118.1 * DEG, 1.25); _append(atoms, oe2)
    elif restype == "HIS":
        cg = cb @ _make_tab(chi[0], 113.6 * DEG, 1.50); _append(atoms, cg)
        nd1 = cg @ _make_tab(chi[1], 122.7 * DEG, 1.38); _append(atoms, nd1)
        cd2 = cg @ _make_tab(chi[1] + 179.9 * DEG, 131.0 * DEG, 1.36); _append(atoms, cd2)
        ce1 = nd1 @ _make_tab(179.9 * DEG, 109.2 * DEG, 1.32); _append(atoms, ce1)
        ne2 = cd2 @ _make_tab(-179.9 * DEG, 107.2 * DEG, 1.37); _append(atoms, ne2)
    elif restype == "ILE":
        cg1 = cb @ _make_tab(chi[0], 110.4 * DEG, 1.53); _append(atoms, cg1)
        cg2 = cb @ _make_tab(chi[0] - 123.2 * DEG, 110.7 * DEG, 1.53); _append(atoms, cg2)
        cd1 = cg1 @ _make_tab(chi[1], 114.0 * DEG, 1.52); _append(atoms, cd1)
    elif restype == "LEU":
        cg = cb @ _make_tab(chi[0], 116.4 * DEG, 1.53); _append(atoms, cg)
        cd1 = cg @ _make_tab(chi[1], 110.4 * DEG, 1.53); _append(atoms, cd1)
        cd2 = cg @ _make_tab(chi[1] + 122.9 * DEG, 110.6 * DEG, 1.53); _append(atoms, cd2)
    elif restype == "LYS":
        cg = cb @ _make_tab(chi[0], 114.0 * DEG, 1.52); _append(atoms, cg)
        cd = cg @ _make_tab(chi[1], 111.5 * DEG, 1.52); _append(atoms, cd)
        ce = cd @ _make_tab(chi[2], 111.6 * DEG, 1.52); _append(atoms, ce)
        nz = ce @ _make_tab(chi[3], 111.8 * DEG, 1.49); _append(atoms, nz)
    elif restype == "MET":
        cg = cb @ _make_tab(chi[0], 113.9 * DEG, 1.52); _append(atoms, cg)
        sd = cg @ _make_tab(chi[1], 112.7 * DEG, 1.81); _append(atoms, sd)
        ce = sd @ _make_tab(chi[2], 100.7 * DEG, 1.79); _append(atoms, ce)
    elif restype == "PHE":
        cg = cb @ _make_tab(chi[0], 113.8 * DEG, 1.50); _append(atoms, cg)
        cd1 = cg @ _make_tab(chi[1], 120.7 * DEG, 1.39); _append(atoms, cd1)
        cd2 = cg @ _make_tab(chi[1] - 180.0 * DEG, 120.5 * DEG, 1.39); _append(atoms, cd2)
        ce1 = cd1 @ _make_tab(-180.0 * DEG, 120.8 * DEG, 1.39); _append(atoms, ce1)
        ce2 = cd2 @ _make_tab(180.0 * DEG, 120.8 * DEG, 1.39); _append(atoms, ce2)
        cz = ce1 @ _make_tab(-0.0 * DEG, 119.9 * DEG, 1.39); _append(atoms, cz)
    elif restype == "PRO":
        cg = cb @ _make_tab(chi[0], 104.2 * DEG, 1.50); _append(atoms, cg)
        cd = cg @ _make_tab(chi[1], 104.9 * DEG, 1.51); _append(atoms, cd)
    elif restype == "SER":
        og = cb @ _make_tab(chi[0], 110.8 * DEG, 1.42); _append(atoms, og)
    elif restype == "THR":
        og1 = cb @ _make_tab(chi[0], 109.2 * DEG, 1.43); _append(atoms, og1)
        cg2 = cb @ _make_tab(chi[0] - 120.4 * DEG, 111.1 * DEG, 1.53); _append(atoms, cg2)
    elif restype == "TRP":
        cg = cb @ _make_tab(chi[0], 113.9 * DEG, 1.50); _append(atoms, cg)
        cd1 = cg @ _make_tab(chi[1], 127.1 * DEG, 1.37); _append(atoms, cd1)
        cd2 = cg @ _make_tab(chi[1] - 179.7 * DEG, 126.6 * DEG, 1.43); _append(atoms, cd2)
        ne1 = cd1 @ _make_tab(-179.8 * DEG, 110.1 * DEG, 1.38); _append(atoms, ne1)
        ce2 = cd2 @ _make_tab(179.8 * DEG, 107.2 * DEG, 1.41); _append(atoms, ce2)
        ce3 = cd2 @ _make_tab(-0.2 * DEG, 133.9 * DEG, 1.40); _append(atoms, ce3)
        cz2 = ce2 @ _make_tab(180.0 * DEG, 122.4 * DEG, 1.40); _append(atoms, cz2)
        cz3 = ce3 @ _make_tab(-180.0 * DEG, 118.7 * DEG, 1.39); _append(atoms, cz3)
        ch2 = cz2 @ _make_tab(-0.0 * DEG, 117.5 * DEG, 1.37); _append(atoms, ch2)
    elif restype == "TYR":
        cg = cb @ _make_tab(chi[0], 113.7 * DEG, 1.51); _append(atoms, cg)
        cd1 = cg @ _make_tab(chi[1], 120.9 * DEG, 1.39); _append(atoms, cd1)
        cd2 = cg @ _make_tab(chi[1] - 179.9 * DEG, 120.8 * DEG, 1.39); _append(atoms, cd2)
        ce1 = cd1 @ _make_tab(-179.9 * DEG, 121.1 * DEG, 1.39); _append(atoms, ce1)
        ce2 = cd2 @ _make_tab(179.9 * DEG, 121.1 * DEG, 1.39); _append(atoms, ce2)
        cz = ce1 @ _make_tab(-0.0 * DEG, 119.5 * DEG, 1.38); _append(atoms, cz)
        oh = cz @ _make_tab(180.0 * DEG, 119.8 * DEG, 1.38); _append(atoms, oh)
    elif restype == "VAL":
        cg1 = cb @ _make_tab(chi[0], 110.7 * DEG, 1.53); _append(atoms, cg1)
        cg2 = cb @ _make_tab(chi[0] + 122.9 * DEG, 110.4 * DEG, 1.53); _append(atoms, cg2)
    else:
        raise KeyError("Unsupported residue type %s" % restype)

    return np.asarray(atoms, dtype="f4")


def build_sidechain_atoms(restype, chi):
    restype = canonical_restype(restype)
    params = RESIDUE_ATOMS[restype]
    if not params:
        return {
            "atom_name": [],
            "atom_element": [],
            "atom_local_pos": np.zeros((0, 3), dtype="f4"),
            "atom_rmin_half": np.zeros(0, dtype="f4"),
            "atom_epsilon": np.zeros(0, dtype="f4"),
            "atom_charge": np.zeros(0, dtype="f4"),
            "bond_index": np.zeros((0, 2), dtype="i4"),
        }

    full_pos = _build_full_positions(restype, chi, psi=0.0)
    side_pos = full_pos[4:]
    if len(params) != side_pos.shape[0]:
        raise RuntimeError("Sidechain geometry/parameter mismatch for %s" % restype)

    return {
        "atom_name": [p["name"] for p in params],
        "atom_element": [p["element"] for p in params],
        "atom_local_pos": side_pos.astype("f4"),
        "atom_rmin_half": np.asarray([p["rmin_half"] for p in params], dtype="f4"),
        "atom_epsilon": np.asarray([p["epsilon"] for p in params], dtype="f4"),
        "atom_charge": np.asarray([p["charge"] for p in params], dtype="f4"),
        "bond_index": np.asarray(RESIDUE_BONDS[restype], dtype="i4"),
    }


def build_sidechain_geometry(restype, chi):
    atom_block = build_sidechain_atoms(restype, chi)
    return {
        "atom_name": atom_block["atom_name"],
        "atom_element": atom_block["atom_element"],
        "atom_local_pos": atom_block["atom_local_pos"],
    }
