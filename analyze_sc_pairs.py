#!/usr/bin/env python3
"""
Compare Upside SC-SC pair interactions against AMBER14SB vacuum + SASA + Env corrections.
"""
import numpy as np
import h5py
import csv
import sys

REPO = "/Users/yinhan/Documents/upside2-md"
E_UP_TO_KCAL = 2.914952774272 / 4.184  # 1 E_up in kcal/mol

# ---------------------------------------------------------------------------
# Residue order (from sidechain.h5 restype_order)
# ---------------------------------------------------------------------------
AA_ORDER = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY',
            'HIS','ILE','LEU','LYS','MET','PHE','PRO','SER',
            'THR','TRP','TYR','VAL']

# ---------------------------------------------------------------------------
# AMBER LJ params: (eps kcal/mol, R_min/2 Å)
# ---------------------------------------------------------------------------
AMBER_LJ = {
    'CT': (0.1094, 1.9080),
    'CA': (0.0860, 1.9080),
    'C':  (0.0860, 1.9080),
    'O':  (0.2100, 1.6612),
    'OH': (0.2104, 1.7210),
    'N':  (0.1700, 1.8240),
    'N2': (0.1700, 1.8240),
    'N3': (0.1700, 1.8240),
    'S':  (0.2500, 2.0000),
    'SH': (0.2500, 2.0000),
}

SASA_RADII = {'CT':1.70,'CA':1.70,'C':1.70,'O':1.52,'OH':1.52,
              'N':1.55,'N2':1.55,'N3':1.55,'S':1.80,'SH':1.80}
PROBE_R = 1.40
GAMMA = 0.0072  # kcal/mol/Å²

# ---------------------------------------------------------------------------
# Geometry builder
# ---------------------------------------------------------------------------
def build_geometry(aa):
    """Return (coords Nx3, types list) for sidechain heavy atoms. COM-centered."""
    cc_sp3 = 1.522
    cc_aro = 1.395
    cn_sp3 = 1.47
    cn_ami = 1.33
    co     = 1.43
    cs     = 1.81
    tet    = np.radians(109.5)
    # Helper: place atom at given distance/angle from previous two
    def extend(prev2, prev1, d, angle, dihedral=180.0):
        b1 = prev1 - prev2
        b1 /= np.linalg.norm(b1)
        # perp in plane
        n_tmp = np.array([1.,0.,0.]) if abs(b1[0]) < 0.9 else np.array([0.,1.,0.])
        perp = np.cross(b1, n_tmp); perp /= np.linalg.norm(perp)
        n2 = np.cross(b1, perp)
        phi = np.radians(dihedral)
        perp2 = np.cos(phi)*perp + np.sin(phi)*n2
        theta = np.pi - angle
        d_vec = np.cos(theta)*b1 + np.sin(theta)*perp2
        return prev1 + d * d_vec

    def branch(base, parent, d, angle, dihedral):
        return extend(parent, base, d, angle, dihedral)

    coords = []
    types  = []

    if aa == 'ALA':
        coords = [[0,0,0]]
        types  = ['CT']
    elif aa == 'ARG':
        # CB-CG-CD-NE-CZ, NH1, NH2
        p0 = np.array([0.,0.,0.])
        p1 = np.array([cc_sp3,0.,0.])
        p2 = extend(p0,p1,cc_sp3,tet)
        p3 = extend(p1,p2,cc_sp3,tet)
        p4 = extend(p2,p3,cn_sp3,tet)   # NE
        p5 = extend(p3,p4,cc_aro,tet)   # CZ (aromatic-ish)
        p6 = extend(p4,p5,cn_ami,tet,120.)  # NH1
        p7 = extend(p4,p5,cn_ami,tet,-120.) # NH2
        coords = [p0,p1,p2,p3,p4,p5,p6,p7]
        types  = ['CT','CT','CT','N2','CA','N2','N2','N2']
        # Only 7 atoms per spec: CB CG CD NE CZ NH1 NH2 = indices 0..6 (drop p7 or include)
        # Spec says 7 atoms [CT,CT,CT,N2,CA,N2,N2] - that's CB,CG,CD,NE,CZ,NH1,NH2
        coords = [p0,p1,p2,p3,p4,p5,p6]
        types  = ['CT','CT','CT','N2','CA','N2','N2']
    elif aa == 'ASN':
        # CB-CG-OD1/ND2
        p0 = np.array([0.,0.,0.])
        p1 = extend(np.array([-cc_sp3,0.,0.]),p0,1.52,tet)  # CG
        p2 = extend(p0,p1,1.23,np.radians(120.),0.)          # OD1 (C=O)
        p3 = extend(p0,p1,cn_ami,np.radians(120.),180.)       # ND2
        coords = [p0,p1,p2,p3]
        types  = ['CT','C','O','N']
    elif aa == 'ASP':
        p0 = np.array([0.,0.,0.])
        p1 = extend(np.array([-cc_sp3,0.,0.]),p0,1.52,tet)
        p2 = extend(p0,p1,1.25,np.radians(120.),0.)
        p3 = extend(p0,p1,1.25,np.radians(120.),180.)
        coords = [p0,p1,p2,p3]
        types  = ['CT','C','O','O']
    elif aa == 'CYS':
        p0 = np.array([0.,0.,0.])
        p1 = extend(np.array([-cc_sp3,0.,0.]),p0,cs,tet)
        coords = [p0,p1]
        types  = ['CT','SH']
    elif aa == 'GLN':
        p0 = np.array([0.,0.,0.])
        p1 = extend(np.array([-cc_sp3,0.,0.]),p0,cc_sp3,tet)
        p2 = extend(p0,p1,1.52,tet)   # CD (carbonyl C)
        p3 = extend(p1,p2,1.23,np.radians(120.),0.)
        p4 = extend(p1,p2,cn_ami,np.radians(120.),180.)
        coords = [p0,p1,p2,p3,p4]
        types  = ['CT','CT','C','O','N']
    elif aa == 'GLU':
        p0 = np.array([0.,0.,0.])
        p1 = extend(np.array([-cc_sp3,0.,0.]),p0,cc_sp3,tet)
        p2 = extend(p0,p1,1.52,tet)
        p3 = extend(p1,p2,1.25,np.radians(120.),0.)
        p4 = extend(p1,p2,1.25,np.radians(120.),180.)
        coords = [p0,p1,p2,p3,p4]
        types  = ['CT','CT','C','O','O']
    elif aa == 'GLY':
        coords = []
        types  = []
    elif aa == 'HIS':
        # CB + imidazole: CG,ND1,CD2,CE1,NE2
        p0 = np.array([0.,0.,0.])  # CB
        p1 = extend(np.array([-cc_sp3,0.,0.]),p0,1.50,tet)  # CG
        # imidazole ring in plane
        p2 = extend(p0,p1,1.38,np.radians(126.),0.)   # ND1
        p3 = extend(p0,p1,1.35,np.radians(126.),180.) # CD2
        p4 = extend(p1,p2,1.32,np.radians(108.),0.)   # CE1
        p5 = extend(p1,p3,1.35,np.radians(108.),0.)   # NE2
        coords = [p0,p1,p2,p3,p4,p5]
        types  = ['CT','CA','N','CA','CA','N']
    elif aa == 'ILE':
        # CB-CG1-CD1 + CG2 branch
        p0 = np.array([0.,0.,0.])
        p1 = extend(np.array([-cc_sp3,0.,0.]),p0,cc_sp3,tet)  # CG1
        p2 = extend(p0,p1,cc_sp3,tet)  # CD1
        p3 = branch(p0,np.array([-cc_sp3,0.,0.]),cc_sp3,tet,120.)  # CG2
        coords = [p0,p1,p2,p3]
        types  = ['CT','CT','CT','CT']
    elif aa == 'LEU':
        p0 = np.array([0.,0.,0.])  # CB
        p1 = extend(np.array([-cc_sp3,0.,0.]),p0,cc_sp3,tet)  # CG
        p2 = extend(p0,p1,cc_sp3,tet,60.)   # CD1
        p3 = extend(p0,p1,cc_sp3,tet,-60.)  # CD2
        coords = [p0,p1,p2,p3]
        types  = ['CT','CT','CT','CT']
    elif aa == 'LYS':
        p0 = np.array([0.,0.,0.])
        p1 = extend(np.array([-cc_sp3,0.,0.]),p0,cc_sp3,tet)
        p2 = extend(p0,p1,cc_sp3,tet)
        p3 = extend(p1,p2,cc_sp3,tet)
        p4 = extend(p2,p3,cn_sp3,tet)  # NZ
        coords = [p0,p1,p2,p3,p4]
        types  = ['CT','CT','CT','CT','N3']
    elif aa == 'MET':
        p0 = np.array([0.,0.,0.])
        p1 = extend(np.array([-cc_sp3,0.,0.]),p0,cs,tet)  # SD
        p2 = extend(p0,p1,cs,tet)  # CE
        coords = [p0,p1,p2]
        types  = ['CT','S','CT']
        # spec: CB-CG-SD-CE = 4 atoms [CT,CT,S,CT]
        p0b = np.array([0.,0.,0.])
        p1b = extend(np.array([-cc_sp3,0.,0.]),p0b,cc_sp3,tet)
        p2b = extend(p0b,p1b,cs,tet)
        p3b = extend(p1b,p2b,cs,tet)
        coords = [p0b,p1b,p2b,p3b]
        types  = ['CT','CT','S','CT']
    elif aa == 'PHE':
        p0 = np.array([0.,0.,0.])  # CB
        p1 = extend(np.array([-cc_sp3,0.,0.]),p0,1.51,tet)  # CG
        # benzene ring: 6 carbons
        ring = []
        for i in range(6):
            angle = np.radians(i * 60.)
            ring.append(p1 + cc_aro * np.array([np.cos(angle+np.pi/2), np.sin(angle+np.pi/2), 0.]))
        coords = [p0, p1] + ring
        types  = ['CT','CA','CA','CA','CA','CA','CA','CA']
        # spec says 7 atoms: CB+CG+CD1+CE1+CZ+CE2+CD2 = CB(CT) + 6 CA
        coords = [p0] + ring
        types  = ['CT','CA','CA','CA','CA','CA','CA']
    elif aa == 'PRO':
        p0 = np.array([0.,0.,0.])
        p1 = extend(np.array([-cc_sp3,0.,0.]),p0,cc_sp3,tet)
        p2 = extend(p0,p1,cc_sp3,tet)
        coords = [p0,p1,p2]
        types  = ['CT','CT','CT']
    elif aa == 'SER':
        p0 = np.array([0.,0.,0.])
        p1 = extend(np.array([-cc_sp3,0.,0.]),p0,co,tet)
        coords = [p0,p1]
        types  = ['CT','OH']
    elif aa == 'THR':
        p0 = np.array([0.,0.,0.])  # CB
        p1 = extend(np.array([-cc_sp3,0.,0.]),p0,co,tet)  # OG1
        p2 = branch(p0,np.array([-cc_sp3,0.,0.]),cc_sp3,tet,120.)  # CG2
        coords = [p0,p1,p2]
        types  = ['CT','OH','CT']
    elif aa == 'TRP':
        # CB + indole: CG,CD1,NE1,CE2,CD2,CE3,CZ3,CH2,CZ2
        p0 = np.array([0.,0.,0.])  # CB
        p1 = extend(np.array([-cc_sp3,0.,0.]),p0,1.50,tet)  # CG
        # Build indole as flat ring system in xy plane
        # 5-ring: CG, CD1, NE1, CE2, CD2
        # 6-ring: CD2, CE3, CZ3, CH2, CZ2, CE2
        # Simplified: place atoms at ring positions
        angle5 = np.radians(108.)
        angle6 = np.radians(120.)
        atoms = [p0, p1]
        # 5-membered ring
        cd1 = p1 + cc_aro*np.array([np.cos(np.pi/5), np.sin(np.pi/5), 0.])
        ne1 = cd1 + cn_ami*np.array([np.cos(3*np.pi/5), np.sin(3*np.pi/5), 0.])
        ce2 = ne1 + cc_aro*np.array([np.cos(np.pi), 0., 0.])
        cd2 = p1 + cc_aro*np.array([np.cos(-np.pi/5), np.sin(-np.pi/5), 0.])
        # 6-membered ring
        ce3 = cd2 + cc_aro*np.array([np.cos(-2*np.pi/6), np.sin(-2*np.pi/6), 0.])
        cz3 = ce3 + cc_aro*np.array([np.cos(-np.pi), 0., 0.])
        ch2 = cz3 + cc_aro*np.array([np.cos(2*np.pi/6), np.sin(2*np.pi/6), 0.])
        cz2 = ce2 + cc_aro*np.array([0., -cc_aro, 0.])
        coords = [p0, p1, cd1, ne1, ce2, cd2, ce3, cz3, ch2, cz2]
        types  = ['CT','CA','CA','N','CA','CA','CA','CA','CA','CA']
    elif aa == 'TYR':
        # CB + phenol ring + OH
        p0 = np.array([0.,0.,0.])  # CB
        ring = []
        cg = p0 + np.array([1.51, 0., 0.])
        for i in range(6):
            angle = np.radians(i * 60.)
            ring.append(cg + cc_aro * np.array([np.cos(angle+np.pi/2), np.sin(angle+np.pi/2), 0.]))
        # OH at para position (index 3)
        oh_pos = ring[3] + np.array([0., 1.36, 0.])
        # spec: CB(CT) + CG,CD1,CE1,CZ,OH,CE2,CD2 -> CT,CA,CA,CA,C,OH,CA,CA = 8 atoms
        coords = [p0] + ring[:3] + [ring[3]] + [oh_pos] + ring[4:]
        types  = ['CT','CA','CA','CA','C','OH','CA','CA']
    elif aa == 'VAL':
        p0 = np.array([0.,0.,0.])  # CB
        p1 = branch(p0,np.array([-cc_sp3,0.,0.]),cc_sp3,tet,60.)   # CG1
        p2 = branch(p0,np.array([-cc_sp3,0.,0.]),cc_sp3,tet,-60.)  # CG2
        coords = [p0,p1,p2]
        types  = ['CT','CT','CT']
    else:
        coords = []
        types  = []

    if len(coords) == 0:
        return np.zeros((0,3)), []

    coords = np.array([np.array(c, dtype=float) for c in coords])
    # COM center
    com = coords.mean(axis=0)
    coords -= com
    return coords, types

# ---------------------------------------------------------------------------
# AMBER LJ energy (vectorized)
# ---------------------------------------------------------------------------
def amber_lj_energy(coords1, types1, coords2, types2):
    """Total LJ energy between two sets of atoms (kcal/mol)."""
    if len(coords1) == 0 or len(coords2) == 0:
        return 0.0
    eps1 = np.array([AMBER_LJ[t][0] for t in types1])
    r1   = np.array([AMBER_LJ[t][1] for t in types1])
    eps2 = np.array([AMBER_LJ[t][0] for t in types2])
    r2   = np.array([AMBER_LJ[t][1] for t in types2])

    # Pairwise distances: (N1, N2, 3) -> (N1, N2)
    diff = coords1[:,None,:] - coords2[None,:,:]  # (N1,N2,3)
    dist = np.linalg.norm(diff, axis=2)            # (N1,N2)

    eps_mix = np.sqrt(eps1[:,None] * eps2[None,:])   # (N1,N2)
    r_mix   = r1[:,None] + r2[None,:]                # (N1,N2) = R_min

    ratio = r_mix / np.maximum(dist, 0.1)
    E = eps_mix * (ratio**12 - 2.*ratio**6)
    return float(np.sum(np.clip(E, -500., 500.)))

# ---------------------------------------------------------------------------
# SASA (Fibonacci sphere)
# ---------------------------------------------------------------------------
def compute_sasa(coords, types):
    """Numerical SASA using Fibonacci sphere probe rolling."""
    if len(coords) == 0:
        return 0.0
    N_pts = 300
    phi = np.pi * (3. - np.sqrt(5.))
    idx = np.arange(N_pts)
    y   = 1. - idx / (N_pts - 1.) * 2.
    r_s = np.sqrt(1. - y*y)
    theta = phi * idx
    pts = np.column_stack([r_s*np.cos(theta), y, r_s*np.sin(theta)])  # (N_pts,3)

    radii = np.array([SASA_RADII[t] + PROBE_R for t in types])  # (N,)
    total = 0.
    for i,(c,rad) in enumerate(zip(coords, radii)):
        sphere_pts = c + rad * pts  # (N_pts,3)
        # Check accessibility: not buried by any other atom
        diffs = sphere_pts[:,None,:] - coords[None,:,:]   # (N_pts,N,3)
        dists = np.linalg.norm(diffs, axis=2)             # (N_pts,N)
        other_rad = radii.copy()
        # mask self
        mask = np.ones(len(coords), dtype=bool)
        mask[i] = False
        buried = np.any(dists[:,mask] < other_rad[None,mask], axis=1)  # (N_pts,)
        n_acc = np.sum(~buried)
        area = 4. * np.pi * (SASA_RADII[types[i]] + PROBE_R)**2 * n_acc / N_pts
        total += area
    return total

def delta_sasa_interaction(coords1, types1, coords2, types2, r_arr):
    """
    Sphere-cap approximation: ΔSASA(r) for bringing mol2 from infinity to distance r.
    We translate mol2 along z-axis by r, compute combined SASA - (SASA1 + SASA2).
    For efficiency: compute SASA at each r value.
    """
    sasa1 = compute_sasa(coords1, types1)
    sasa2 = compute_sasa(coords2, types2)
    sasa_inf = sasa1 + sasa2
    result = []
    for r in r_arr:
        # Orient mol2 along +z at distance r
        c2_shift = coords2 + np.array([0., 0., r])
        combined_coords = np.vstack([coords1, c2_shift])
        combined_types = types1 + types2
        sasa_r = compute_sasa(combined_coords, combined_types)
        result.append(GAMMA * (sasa_r - sasa_inf))
    return np.array(result)

# ---------------------------------------------------------------------------
# Rotation-averaged LJ energy
# ---------------------------------------------------------------------------
def random_rotation_matrices(N, rng):
    """Generate N random rotation matrices via QR decomposition."""
    A = rng.standard_normal((N, 3, 3))
    Q, R = np.linalg.qr(A)
    # Fix sign to ensure det=+1
    signs = np.sign(np.linalg.det(Q))
    Q *= signs[:,None,None]
    return Q.astype(np.float32)

def rot_avg_lj(coords1, types1, coords2, types2, r_arr, N_ROT=600):
    """Rotation-averaged LJ energy vs r (kcal/mol)."""
    if len(coords1) == 0 or len(coords2) == 0:
        return np.zeros(len(r_arr))

    rng = np.random.default_rng(42)
    Rmats = random_rotation_matrices(N_ROT, rng)  # (N_ROT,3,3)

    eps1 = np.array([AMBER_LJ[t][0] for t in types1], dtype=np.float32)
    r1v  = np.array([AMBER_LJ[t][1] for t in types1], dtype=np.float32)
    eps2 = np.array([AMBER_LJ[t][0] for t in types2], dtype=np.float32)
    r2v  = np.array([AMBER_LJ[t][1] for t in types2], dtype=np.float32)

    eps_mix = np.sqrt(eps1[:,None] * eps2[None,:])  # (N1,N2)
    r_mix   = r1v[:,None] + r2v[None,:]              # (N1,N2)

    c1f = coords1.astype(np.float32)  # (N1,3)
    c2f = coords2.astype(np.float32)  # (N2,3)

    results = []
    for r in r_arr:
        # Rotate mol1 (fixed) and mol2 (offset by r along z after rotation)
        # Rotate mol2 with all N_ROT rotations
        c2_rot = np.einsum('rij,nj->rni', Rmats, c2f)  # (N_ROT,N2,3)
        c2_rot[:,:,2] += r   # shift along z

        # Distances: (N_ROT, N1, N2)
        # Use broadcasting: c1f (N1,3), c2_rot (N_ROT,N2,3)
        diff = c1f[None,:,None,:] - c2_rot[:,None,:,:]  # (N_ROT,N1,N2,3)
        dist = np.sqrt(np.sum(diff**2, axis=3))          # (N_ROT,N1,N2)

        ratio = r_mix[None,:,:] / np.maximum(dist, 0.1)
        E = eps_mix[None,:,:] * (ratio**12 - 2.*ratio**6)
        E = np.clip(E, -500., 500.)
        E_per_rot = E.sum(axis=(1,2))  # (N_ROT,)
        results.append(float(E_per_rot.mean()))

    return np.array(results)

# ---------------------------------------------------------------------------
# Upside spline evaluator
# ---------------------------------------------------------------------------
def rr(v): return np.roll(v, 1)

def deBoor(k, x):
    xb = int(x)
    y = x - xb
    c = np.array(k[xb-1:xb+3], dtype=float)
    yu = y - np.array([0., -2., -1., 0.])
    c1 = (1-yu/3)*rr(c) + (yu/3)*c
    c2 = (1-yu/2)*rr(c1) + (yu/2)*c1
    c3 = (1-yu)*rr(c2) + yu*c2
    return float(c3[3])

def cDB(k, x, n):
    if x <= 1.: return deBoor(k, 1.0)
    if x >= n-2: return deBoor(k, float(n-2)-1e-6)
    return deBoor(k, x)

def eval_ang(k, c):
    INV_DTHETA = (15-3)/2.  # = 6.0
    x = np.clip((c+1)*INV_DTHETA + 1., 1., 15.-2.-1e-6)
    return deBoor(k, x)

def eval_radial_spline(k, r, r0=0., dr=0.7, n=12):
    """Evaluate radial spline at distance r (Å). r0=0, dr=0.7, n=12 knots."""
    x = (r - r0) / dr + 1.
    return cDB(k, x, n)

# Angular average of angular spline: mean over cos_theta uniform in [-1,1]
N_ANG = 200
COS_VALS = np.linspace(-1., 1., N_ANG)

def ang_avg(k):
    return np.mean([eval_ang(k, c) for c in COS_VALS])

def upside_pair_energy(knots, r_arr):
    """Evaluate orientation-averaged Upside pair energy (E_up) vs r."""
    ang1_k = knots[0:15]
    ang2_k = knots[15:30]
    wide_k = knots[30:42]
    narr_k = knots[42:54]

    a1 = ang_avg(ang1_k)
    a2 = ang_avg(ang2_k)

    energies = []
    for r in r_arr:
        wide = eval_radial_spline(wide_k, r)
        narr = eval_radial_spline(narr_k, r)
        e = wide + a1 * a2 * narr
        energies.append(e)
    return np.array(energies) * E_UP_TO_KCAL

# ---------------------------------------------------------------------------
# Environment correction
# ---------------------------------------------------------------------------
def compact_sigmoid(x, s):
    z = x * s
    z = np.asarray(z, dtype=float)
    result = np.where(z <= -1., 0., np.where(z >= 1., 1., 0.5 + 0.5*z*(1.5 - z*z*0.5)))
    return result

# Pre-compute angular average of compact_sigmoid(0.3 - cos_theta, 0.6) over cos_theta in [-1,1]
_cos_grid = np.linspace(-1., 1., 500)
ANG_AVG_CS = float(np.mean(compact_sigmoid(0.3 - _cos_grid, 0.6)))

def env_cov(r):
    """Coverage of bead i from bead j at distance r (scalar or array)."""
    radial = compact_sigmoid(r - 8.0, 1.0)
    return radial * ANG_AVG_CS

def env_delta(r_arr, type_idx, scale, center, sharpness):
    """ΔE_env (E_up) for residue type_idx as function of r."""
    cov_r = env_cov(np.asarray(r_arr, dtype=float))
    cov_0 = 0.  # at r=inf, coverage = 0
    dE = scale[type_idx] * (
        compact_sigmoid(center[type_idx] - cov_r, sharpness[type_idx]) -
        compact_sigmoid(center[type_idx] - cov_0, sharpness[type_idx])
    )
    return dE  # E_up

def env_two_body(r_arr, ti, tj, scale, center, sharpness):
    """Total two-body env correction in kcal/mol."""
    dEi = env_delta(r_arr, ti, scale, center, sharpness)
    dEj = env_delta(r_arr, tj, scale, center, sharpness)
    return (dEi + dEj) * E_UP_TO_KCAL

# ---------------------------------------------------------------------------
# Load parameters
# ---------------------------------------------------------------------------
print("Loading parameters...")
with h5py.File(f"{REPO}/parameters/ff_2.1/sidechain.h5", 'r') as f:
    pair_interaction = f['pair_interaction'][:]  # (20,20,54)
    restype_order    = [x.decode() for x in f['restype_order'][:]]

with h5py.File(f"{REPO}/parameters/ff_2.1/environment.h5", 'r') as f:
    env_scale     = f['scale'][:]     # (20,)
    env_center    = f['center'][:]    # (20,)
    env_sharpness = f['sharpness'][:] # (20,)

# Map AA_ORDER index to restype_order index
aa_to_rt = {aa: restype_order.index(aa) for aa in AA_ORDER}

print("Building geometries...")
geometries = {}
for aa in AA_ORDER:
    coords, types = build_geometry(aa)
    geometries[aa] = (coords, types)

# r grid for evaluation
r_arr = np.arange(2.0, 14.1, 0.2)

print(f"Evaluating 210 pairs on {len(r_arr)} r points...")

def find_well(r_arr, E_arr):
    """Return (E_min, r_at_min)."""
    if len(E_arr) == 0 or np.all(E_arr == 0):
        return 0., float('nan')
    idx = np.argmin(E_arr)
    return float(E_arr[idx]), float(r_arr[idx])

rows = []
pair_count = 0
total_pairs = sum(1 for i in range(20) for j in range(i, 20))

for i, aa1 in enumerate(AA_ORDER):
    for j in range(i, 20):
        aa2 = AA_ORDER[j]
        pair_count += 1

        c1, t1 = geometries[aa1]
        c2, t2 = geometries[aa2]

        ri = aa_to_rt[aa1]
        rj = aa_to_rt[aa2]

        # Upside pair potential
        knots = pair_interaction[ri, rj]
        E_up = upside_pair_energy(knots, r_arr)
        e_up, r_up = find_well(r_arr, E_up)

        # GLY has no sidechain - AMBER vac = 0
        if len(t1) == 0 or len(t2) == 0:
            e_amb, r_amb = 0., float('nan')
            e_sa,  r_sa  = 0., float('nan')
            e_env, r_env = 0., float('nan')
            rows.append({
                'Pair': f"{aa1}-{aa2}",
                'AMBvac': e_amb, 'rAMB': r_amb,
                'AMB+SA': e_sa,  'rSA':  r_sa,
                'AMB+SA+Env': e_env, 'rEnv': r_env,
                'Upside': e_up,  'rU':   r_up,
            })
            if pair_count % 20 == 0:
                print(f"  {pair_count}/{total_pairs} pairs done")
            continue

        # AMBER vacuum (rotation-averaged)
        E_amb = rot_avg_lj(c1, t1, c2, t2, r_arr, N_ROT=600)
        e_amb, r_amb = find_well(r_arr, E_amb)

        # SASA correction (slower, compute at subset of r)
        r_sa_grid = np.arange(2.0, 14.1, 0.5)
        dSASA = delta_sasa_interaction(c1, t1, c2, t2, r_sa_grid)
        # Interpolate to full r_arr
        dSASA_full = np.interp(r_arr, r_sa_grid, dSASA)
        E_amb_sa = E_amb + dSASA_full
        e_sa, r_sa = find_well(r_arr, E_amb_sa)

        # Environment two-body correction
        dE_env = env_two_body(r_arr, ri, rj, env_scale, env_center, env_sharpness)
        E_amb_sa_env = E_amb_sa + dE_env
        e_env, r_env = find_well(r_arr, E_amb_sa_env)

        rows.append({
            'Pair': f"{aa1}-{aa2}",
            'AMBvac': e_amb, 'rAMB': r_amb,
            'AMB+SA': e_sa,  'rSA':  r_sa,
            'AMB+SA+Env': e_env, 'rEnv': r_env,
            'Upside': e_up,  'rU':   r_up,
        })

        if pair_count % 10 == 0:
            print(f"  {pair_count}/{total_pairs} pairs done")

# ---------------------------------------------------------------------------
# Print table
# ---------------------------------------------------------------------------
print()
hdr = f"{'Pair':<10} {'AMBvac':>10} {'rAMB':>7} {'AMB+SA':>10} {'rSA':>7} {'AMB+SA+Env':>12} {'rEnv':>7} {'Upside':>10} {'rU':>7}"
print(hdr)
print('-'*len(hdr))
for row in rows:
    def fmt(v): return f"{v:7.3f}" if not (isinstance(v,float) and np.isnan(v)) else "    NaN"
    print(f"{row['Pair']:<10} {row['AMBvac']:>10.3f} {fmt(row['rAMB'])} {row['AMB+SA']:>10.3f} "
          f"{fmt(row['rSA'])} {row['AMB+SA+Env']:>12.3f} {fmt(row['rEnv'])} "
          f"{row['Upside']:>10.3f} {fmt(row['rU'])}")

# Summary stats (excluding GLY pairs)
non_gly = [r for r in rows if 'GLY' not in r['Pair']]
mad_amb_up = np.mean([abs(r['AMBvac'] - r['Upside']) for r in non_gly])
mad_env_up = np.mean([abs(r['AMB+SA+Env'] - r['Upside']) for r in non_gly])
print()
print(f"Summary (excluding GLY, N={len(non_gly)} pairs):")
print(f"  MAD(AMBvac vs Upside):       {mad_amb_up:.3f} kcal/mol")
print(f"  MAD(AMB+SA+Env vs Upside):   {mad_env_up:.3f} kcal/mol")

# Save CSV
csv_path = "/tmp/sc_pair_comparison.csv"
fieldnames = ['Pair','AMBvac','rAMB','AMB+SA','rSA','AMB+SA+Env','rEnv','Upside','rU']
with open(csv_path, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)
print(f"\nResults saved to {csv_path}")
