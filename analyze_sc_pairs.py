#!/usr/bin/env python3
"""
Compare Upside SC-SC pair interactions against AMBER14SB and CHARMM36.

PMF: W(r) = -kT ln<exp(-beta*U(r,Omega))>_Omega
This is the proper potential of mean force for coarse-graining a rigid SC group
to a single bead, averaging over all relative orientations at temperature T.
"""
import numpy as np
from scipy.special import logsumexp
import h5py
import csv

REPO = "/Users/yinhan/Documents/upside2-md"
E_UP_TO_KCAL = 2.914952774272 / 4.184  # 1 E_up in kcal/mol
T_REF = 300.                             # K for PMF
kT = 0.001987 * T_REF                   # kcal/mol
beta = 1.0 / kT

AA_ORDER = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY',
            'HIS','ILE','LEU','LYS','MET','PHE','PRO','SER',
            'THR','TRP','TYR','VAL']

# AMBER14SB LJ: (eps kcal/mol, Rmin/2 Å)
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

# CHARMM36 LJ: (eps kcal/mol, Rmin/2 Å)
CHARMM_LJ = {
    'CT1': (0.020, 2.275),
    'CT2': (0.056, 2.175),
    'CT3': (0.080, 2.060),
    'CA':  (0.070, 2.000),
    'CC':  (0.070, 2.000),   # carboxylate/amide carbonyl C
    'CPH1':(0.052, 1.800),   # His ring C (CG, CD2)
    'CPH2':(0.052, 1.800),   # His ring CH (CE1)
    'CW':  (0.076, 1.900),   # Trp 5-ring C
    'CY':  (0.076, 1.900),   # Trp 6-ring C fused
    'NC2': (0.200, 1.850),   # Arg guanidinium N
    'NH2': (0.200, 1.850),   # amide NH2 (Asn, Gln)
    'NH3': (0.200, 1.850),   # Lys NZ
    'NR1': (0.200, 1.850),   # His ND1 (NH)
    'NR2': (0.200, 1.850),   # His NE2 (N:)
    'NY':  (0.200, 1.850),   # Trp NE1
    'O':   (0.120, 1.700),   # amide/carbonyl O
    'OC':  (0.120, 1.700),   # carboxylate O
    'OH1': (0.152, 1.770),   # hydroxyl O
    'S':   (0.450, 2.000),   # CYS SG, MET SD
}

# Bondi radii for SASA (element-based, independent of FF)
SASA_RADII = {
    # AMBER types
    'CT':1.70,'CA':1.70,'C':1.70,'O':1.52,'OH':1.52,
    'N':1.55,'N2':1.55,'N3':1.55,'S':1.80,'SH':1.80,
    # CHARMM types
    'CT1':1.70,'CT2':1.70,'CT3':1.70,'CC':1.70,
    'CPH1':1.70,'CPH2':1.70,'CW':1.70,'CY':1.70,
    'NC2':1.55,'NH2':1.55,'NH3':1.55,'NR1':1.55,'NR2':1.55,'NY':1.55,
    'OC':1.52,'OH1':1.52,
}
PROBE_R = 1.40
GAMMA   = 0.0072  # kcal/mol/Å²

# ---------------------------------------------------------------------------
# Geometry builder  (coords, amber_types, charmm_types)
# ---------------------------------------------------------------------------
def build_geometry(aa):
    """COM-centered heavy-atom sidechain geometry with AMBER and CHARMM36 types."""
    cc = 1.522   # sp3 C-C
    ca = 1.395   # aromatic C-C
    cn = 1.47    # sp3 C-N
    co = 1.43    # C-O
    cs = 1.81    # C-S
    tet = np.radians(109.5)

    def extend(prev2, prev1, d, angle, dihedral=180.0):
        b = prev1 - prev2; b /= np.linalg.norm(b)
        n = np.array([1.,0.,0.]) if abs(b[0]) < 0.9 else np.array([0.,1.,0.])
        p = np.cross(b, n); p /= np.linalg.norm(p)
        n2 = np.cross(b, p)
        phi = np.radians(dihedral)
        perp2 = np.cos(phi)*p + np.sin(phi)*n2
        return prev1 + d*(np.cos(np.pi - angle)*b + np.sin(np.pi - angle)*perp2)

    origin = np.zeros(3)
    back   = np.array([-cc, 0., 0.])

    coords = []
    at     = []   # AMBER types
    ct     = []   # CHARMM types

    if aa == 'ALA':
        coords, at, ct = [origin], ['CT'], ['CT3']

    elif aa == 'ARG':
        p0 = origin
        p1 = np.array([cc,0.,0.])
        p2 = extend(p0,p1,cc,tet)
        p3 = extend(p1,p2,cc,tet)           # NE
        p4 = extend(p2,p3,1.45,tet)         # CZ
        p5 = extend(p3,p4,1.33,np.radians(120.), 120.)  # NH1
        p6 = extend(p3,p4,1.33,np.radians(120.),-120.)  # NH2
        coords = [p0,p1,p2,p3,p4,p5,p6]
        at     = ['CT','CT','CT','N2','CA','N2','N2']
        ct     = ['CT2','CT2','CT2','NC2','CC','NC2','NC2']

    elif aa == 'ASN':
        p0 = origin
        p1 = extend(back,p0,1.52,tet)
        p2 = extend(p0,p1,1.23,np.radians(120.), 0.)
        p3 = extend(p0,p1,1.33,np.radians(120.),180.)
        coords = [p0,p1,p2,p3]
        at     = ['CT','C','O','N']
        ct     = ['CT2','CC','O','NH2']

    elif aa == 'ASP':
        p0 = origin
        p1 = extend(back,p0,1.52,tet)
        p2 = extend(p0,p1,1.25,np.radians(120.), 0.)
        p3 = extend(p0,p1,1.25,np.radians(120.),180.)
        coords = [p0,p1,p2,p3]
        at     = ['CT','C','O','O']
        ct     = ['CT2','CC','OC','OC']

    elif aa == 'CYS':
        p0 = origin
        p1 = extend(back,p0,cs,tet)
        coords = [p0,p1]
        at     = ['CT','SH']
        ct     = ['CT2','S']

    elif aa == 'GLN':
        p0 = origin
        p1 = extend(back,p0,cc,tet)
        p2 = extend(p0,p1,1.52,tet)
        p3 = extend(p1,p2,1.23,np.radians(120.), 0.)
        p4 = extend(p1,p2,1.33,np.radians(120.),180.)
        coords = [p0,p1,p2,p3,p4]
        at     = ['CT','CT','C','O','N']
        ct     = ['CT2','CT2','CC','O','NH2']

    elif aa == 'GLU':
        p0 = origin
        p1 = extend(back,p0,cc,tet)
        p2 = extend(p0,p1,1.52,tet)
        p3 = extend(p1,p2,1.25,np.radians(120.), 0.)
        p4 = extend(p1,p2,1.25,np.radians(120.),180.)
        coords = [p0,p1,p2,p3,p4]
        at     = ['CT','CT','C','O','O']
        ct     = ['CT2','CT2','CC','OC','OC']

    elif aa == 'GLY':
        return np.zeros((0,3)), [], []

    elif aa == 'HIS':
        p0 = origin  # CB
        p1 = extend(back,p0,1.50,tet)  # CG
        p2 = extend(p0,p1,1.38,np.radians(126.), 0.)   # ND1
        p3 = extend(p0,p1,1.35,np.radians(126.),180.)  # CD2
        p4 = extend(p1,p2,1.32,np.radians(108.), 0.)   # CE1
        p5 = extend(p1,p3,1.35,np.radians(108.), 0.)   # NE2
        coords = [p0,p1,p2,p3,p4,p5]
        at     = ['CT','CA','N','CA','CA','N']
        ct     = ['CT2','CPH1','NR1','CPH1','CPH2','NR2']

    elif aa == 'ILE':
        p0 = origin             # CB
        p1 = extend(back,p0,cc,tet)     # CG1
        p2 = extend(p0,p1,cc,tet)       # CD1
        p3 = extend(back,p0,cc,tet,120.)# CG2 (branch at CA side)
        coords = [p0,p1,p2,p3]
        at     = ['CT','CT','CT','CT']
        ct     = ['CT1','CT2','CT3','CT3']

    elif aa == 'LEU':
        p0 = origin             # CB
        p1 = extend(back,p0,cc,tet)    # CG
        p2 = extend(p0,p1,cc,tet, 60.) # CD1
        p3 = extend(p0,p1,cc,tet,-60.) # CD2
        coords = [p0,p1,p2,p3]
        at     = ['CT','CT','CT','CT']
        ct     = ['CT2','CT1','CT3','CT3']

    elif aa == 'LYS':
        p0 = origin
        p1 = extend(back,p0,cc,tet)
        p2 = extend(p0,p1,cc,tet)
        p3 = extend(p1,p2,cc,tet)
        p4 = extend(p2,p3,cn,tet)   # NZ
        coords = [p0,p1,p2,p3,p4]
        at     = ['CT','CT','CT','CT','N3']
        ct     = ['CT2','CT2','CT2','CT2','NH3']

    elif aa == 'MET':
        p0 = origin
        p1 = extend(back,p0,cc,tet)
        p2 = extend(p0,p1,cs,tet)    # SD
        p3 = extend(p1,p2,cs,tet)    # CE
        coords = [p0,p1,p2,p3]
        at     = ['CT','CT','S','CT']
        ct     = ['CT2','CT2','S','CT3']

    elif aa == 'PHE':
        p0 = origin  # CB
        cg = np.array([1.51, 0., 0.])
        ring = [cg + ca*np.array([np.cos(a+np.pi/2), np.sin(a+np.pi/2), 0.])
                for a in [i*np.pi/3 for i in range(6)]]
        coords = [p0] + ring
        at     = ['CT'] + ['CA']*6
        ct     = ['CT2'] + ['CA']*6

    elif aa == 'PRO':
        p0 = origin
        p1 = extend(back,p0,cc,tet)
        p2 = extend(p0,p1,cc,tet)
        coords = [p0,p1,p2]
        at     = ['CT','CT','CT']
        ct     = ['CT2','CT2','CT2']

    elif aa == 'SER':
        p0 = origin
        p1 = extend(back,p0,co,tet)
        coords = [p0,p1]
        at     = ['CT','OH']
        ct     = ['CT2','OH1']

    elif aa == 'THR':
        p0 = origin                         # CB
        p1 = extend(back,p0,co,tet)         # OG1
        p2 = extend(back,p0,cc,tet,120.)    # CG2
        coords = [p0,p1,p2]
        at     = ['CT','OH','CT']
        ct     = ['CT1','OH1','CT3']

    elif aa == 'TRP':
        # CB + indole (5+6 fused ring)
        p0 = origin  # CB
        cg = np.array([1.50, 0., 0.])  # CG
        # 5-membered ring: CG, CD1, NE1, CE2, CD2
        cd1 = cg + 1.37*np.array([np.cos(np.pi/5),  np.sin(np.pi/5),  0.])
        ne1 = cd1 + 1.37*np.array([np.cos(3*np.pi/5), np.sin(3*np.pi/5), 0.])
        ce2 = ne1 + 1.37*np.array([np.cos(np.pi),    0.,               0.])
        cd2 = cg  + 1.44*np.array([np.cos(-np.pi/5), np.sin(-np.pi/5), 0.])
        # 6-membered ring: CD2, CE3, CZ3, CH2, CZ2, CE2
        ce3 = cd2 + ca*np.array([np.cos(-2*np.pi/6), np.sin(-2*np.pi/6), 0.])
        cz3 = ce3 + ca*np.array([np.cos(-np.pi), 0., 0.])
        ch2 = cz3 + ca*np.array([np.cos(2*np.pi/6),  np.sin(2*np.pi/6), 0.])
        cz2 = ce2 + np.array([0., -ca, 0.])
        coords = [p0, cg, cd1, ne1, ce2, cd2, ce3, cz3, ch2, cz2]
        at     = ['CT','CA','CA','N','CA','CA','CA','CA','CA','CA']
        ct     = ['CT2','CW','CW','NY','CY','CY','CA','CA','CA','CY']

    elif aa == 'TYR':
        p0 = origin  # CB
        cg = np.array([1.51, 0., 0.])
        ring = [cg + ca*np.array([np.cos(a+np.pi/2), np.sin(a+np.pi/2), 0.])
                for a in [i*np.pi/3 for i in range(6)]]
        oh = ring[3] + np.array([0., 1.36, 0.])   # OH at para
        coords = [p0] + ring[:4] + [oh] + ring[4:]
        at     = ['CT','CA','CA','CA','CA','OH','CA','CA']
        ct     = ['CT2','CA','CA','CA','CA','OH1','CA','CA']

    elif aa == 'VAL':
        p0 = origin                         # CB
        p1 = extend(back,p0,cc,tet, 60.)    # CG1
        p2 = extend(back,p0,cc,tet,-60.)    # CG2
        coords = [p0,p1,p2]
        at     = ['CT','CT','CT']
        ct     = ['CT1','CT3','CT3']

    else:
        return np.zeros((0,3)), [], []

    coords = np.array([np.array(c, dtype=float) for c in coords])
    coords -= coords.mean(axis=0)  # COM center
    return coords, at, ct


# ---------------------------------------------------------------------------
# PMF  W(r) = -kT ln<exp(-beta*U(r,Omega))>_Omega
# ---------------------------------------------------------------------------
def random_rotation_matrices(N, rng):
    A = rng.standard_normal((N,3,3))
    Q, R = np.linalg.qr(A)
    Q *= np.sign(np.linalg.det(Q))[:,None,None]
    return Q.astype(np.float32)


def pmf_vs_r(coords1, types1, coords2, types2, r_arr, ff_dict, N_ROT=800):
    """
    W(r) = -kT ln[(1/N) sum_i exp(-beta * U(r, Omega_i))]
    Returns array of shape (len(r_arr),) in kcal/mol.
    """
    if len(coords1) == 0 or len(coords2) == 0:
        return np.zeros(len(r_arr))

    rng = np.random.default_rng(42)
    Rmats = random_rotation_matrices(N_ROT, rng)  # (N_ROT,3,3)

    eps1 = np.array([ff_dict[t][0] for t in types1], dtype=np.float32)
    rm1  = np.array([ff_dict[t][1] for t in types1], dtype=np.float32)
    eps2 = np.array([ff_dict[t][0] for t in types2], dtype=np.float32)
    rm2  = np.array([ff_dict[t][1] for t in types2], dtype=np.float32)

    eps_mix = np.sqrt(eps1[:,None] * eps2[None,:])   # (N1,N2)
    r_mix   = rm1[:,None] + rm2[None,:]              # (N1,N2) = R_min

    c1f = coords1.astype(np.float32)
    c2f = coords2.astype(np.float32)

    W = np.empty(len(r_arr))
    for ki, r in enumerate(r_arr):
        # Rotate mol2 with all rotations, then shift along z
        c2_rot = np.einsum('rij,nj->rni', Rmats, c2f)  # (N_ROT,N2,3)
        c2_rot[:,:,2] += r

        # Pairwise distances (N_ROT, N1, N2)
        diff = c1f[None,:,None,:] - c2_rot[:,None,:,:]
        dist = np.sqrt((diff**2).sum(-1))

        # LJ energy per rotation (N_ROT,)
        x   = r_mix[None,:,:] / np.maximum(dist, 0.3)
        E   = (eps_mix[None,:,:] * (x**12 - 2.*x**6)).sum(axis=(1,2))

        # Boltzmann-weighted PMF via log-sum-exp
        logZ = float(logsumexp(-beta * E.astype(np.float64))) - np.log(N_ROT)
        W[ki] = -kT * logZ

    return W


# ---------------------------------------------------------------------------
# SASA (Fibonacci sphere)
# ---------------------------------------------------------------------------
def compute_sasa(coords, types):
    if len(coords) == 0:
        return 0.0
    N_pts = 300
    phi   = np.pi * (3. - np.sqrt(5.))
    idx   = np.arange(N_pts)
    y_s   = 1. - idx / (N_pts - 1.) * 2.
    r_s   = np.sqrt(1. - y_s**2)
    pts   = np.column_stack([r_s*np.cos(phi*idx), y_s, r_s*np.sin(phi*idx)])
    radii = np.array([SASA_RADII[t] + PROBE_R for t in types])
    total = 0.
    for i, (c, rad) in enumerate(zip(coords, radii)):
        sp = c + rad * pts
        d  = np.linalg.norm(sp[:,None,:] - coords[None,:,:], axis=2)
        mask = np.ones(len(coords), dtype=bool); mask[i] = False
        n_acc = np.sum(~np.any(d[:,mask] < radii[None,mask], axis=1))
        total += 4.*np.pi*(SASA_RADII[types[i]]+PROBE_R)**2 * n_acc/N_pts
    return total


def delta_sasa(coords1, types1, coords2, types2, r_arr):
    sasa_inf = compute_sasa(coords1, types1) + compute_sasa(coords2, types2)
    out = []
    for r in r_arr:
        c2s = coords2 + np.array([0., 0., r])
        sasa_r = compute_sasa(np.vstack([coords1, c2s]), types1 + types2)
        out.append(GAMMA * (sasa_r - sasa_inf))
    return np.array(out)


# ---------------------------------------------------------------------------
# Upside spline evaluator (orientation-averaged energy)
# ---------------------------------------------------------------------------
def rr(v): return np.roll(v, 1)

def deBoor(k, x):
    xb = int(x); y = x - xb
    c  = np.array(k[xb-1:xb+3], dtype=float)
    yu = y - np.array([0.,-2.,-1.,0.])
    c1 = (1-yu/3)*rr(c) + (yu/3)*c
    c2 = (1-yu/2)*rr(c1) + (yu/2)*c1
    c3 = (1-yu)*rr(c2) + yu*c2
    return float(c3[3])

def cDB(k, x, n):
    if x <= 1.:    return deBoor(k, 1.0)
    if x >= n-2:   return deBoor(k, float(n-2)-1e-6)
    return deBoor(k, x)

def eval_ang(k, c):
    return cDB(k, (c+1)*6.0+1., 15)

def eval_rad(k, r, dr=0.7, n=12):
    return cDB(k, r/dr + 1., n)

_cos = np.linspace(-1., 1., 300)

def ang_avg(k):
    return float(np.mean([eval_ang(k, c) for c in _cos]))

def upside_pmf(knots, r_arr):
    """Orientation-averaged Upside pair energy (kcal/mol)."""
    a1 = ang_avg(knots[0:15])
    a2 = ang_avg(knots[15:30])
    return np.array([eval_rad(knots[30:42], r) + a1*a2*eval_rad(knots[42:54], r)
                     for r in r_arr]) * E_UP_TO_KCAL


# ---------------------------------------------------------------------------
# Environment two-body correction (tiny, included for completeness)
# ---------------------------------------------------------------------------
def compact_sigmoid(x, s):
    z = np.asarray(x)*s
    return np.where(z<=-1., 0., np.where(z>=1., 1., 0.5+0.5*z*(1.5-z*z*0.5)))

_cg = np.linspace(-1.,1.,500)
_ANG_ENV = float(np.mean(compact_sigmoid(0.3-_cg, 0.6)))

def env_two_body(r_arr, ti, tj, scale, center, sharpness):
    cov = compact_sigmoid(np.asarray(r_arr)-8., 1.) * _ANG_ENV
    dEi = scale[ti]*(compact_sigmoid(center[ti]-cov, sharpness[ti])
                    -compact_sigmoid(center[ti],     sharpness[ti]))
    dEj = scale[tj]*(compact_sigmoid(center[tj]-cov, sharpness[tj])
                    -compact_sigmoid(center[tj],     sharpness[tj]))
    return (dEi + dEj) * E_UP_TO_KCAL


# ---------------------------------------------------------------------------
# Load parameters
# ---------------------------------------------------------------------------
print("Loading parameters...")
with h5py.File(f"{REPO}/parameters/ff_2.1/sidechain.h5", 'r') as f:
    pair_interaction = f['pair_interaction'][:]
    restype_order    = [x.decode() for x in f['restype_order'][:]]

with h5py.File(f"{REPO}/parameters/ff_2.1/environment.h5", 'r') as f:
    env_scale     = f['scale'][:]
    env_center    = f['center'][:]
    env_sharpness = f['sharpness'][:]

aa_to_rt = {aa: restype_order.index(aa) for aa in AA_ORDER}

print("Building geometries...")
geometries = {aa: build_geometry(aa) for aa in AA_ORDER}

r_arr    = np.arange(2.0, 14.1, 0.2)
r_sa_grid = np.arange(2.0, 14.1, 0.5)

print(f"Evaluating 210 pairs on {len(r_arr)} r points (T={T_REF}K PMF)...")

def well(r, E):
    if len(E)==0 or np.all(E==0): return 0., float('nan')
    i = int(np.argmin(E))
    return float(E[i]), float(r[i])

rows = []
for pi, aa1 in enumerate(AA_ORDER):
    for pj in range(pi, 20):
        aa2 = AA_ORDER[pj]
        c1, at1, ct1 = geometries[aa1]
        c2, at2, ct2 = geometries[aa2]
        ri = aa_to_rt[aa1]
        rj = aa_to_rt[aa2]

        # Upside orientation-averaged pair energy
        E_up = upside_pmf(pair_interaction[ri, rj], r_arr)
        e_up, r_up = well(r_arr, E_up)

        if len(at1)==0 or len(at2)==0:
            rows.append({'Pair':f"{aa1}-{aa2}",
                         'AMBvac':0.,'rAMB':float('nan'),
                         'AMB+SA':0.,'rSA': float('nan'),
                         'CHARvac':0.,'rCHAR':float('nan'),
                         'CHAR+SA':0.,'rCSA': float('nan'),
                         'Upside':e_up,'rU':r_up})
            continue

        # AMBER PMF
        W_amb = pmf_vs_r(c1, at1, c2, at2, r_arr, AMBER_LJ)
        e_amb, r_amb = well(r_arr, W_amb)

        # CHARMM36 PMF
        W_chm = pmf_vs_r(c1, ct1, c2, ct2, r_arr, CHARMM_LJ)
        e_chm, r_chm = well(r_arr, W_chm)

        # SASA correction (geometry/element independent of FF)
        dSA = np.interp(r_arr, r_sa_grid,
                        delta_sasa(c1, at1, c2, at2, r_sa_grid))

        # Env two-body correction
        dEnv = env_two_body(r_arr, ri, rj, env_scale, env_center, env_sharpness)

        e_ambsa, r_ambsa   = well(r_arr, W_amb + dSA + dEnv)
        e_chmsa, r_chmsa   = well(r_arr, W_chm + dSA + dEnv)

        rows.append({'Pair':f"{aa1}-{aa2}",
                     'AMBvac': e_amb,   'rAMB':  r_amb,
                     'AMB+SA': e_ambsa, 'rSA':   r_ambsa,
                     'CHARvac':e_chm,   'rCHAR': r_chm,
                     'CHAR+SA':e_chmsa, 'rCSA':  r_chmsa,
                     'Upside': e_up,    'rU':    r_up})

        n_done = sum(1 for r in rows if 'GLY' not in r['Pair'])
        if n_done % 15 == 0:
            print(f"  {len(rows)}/210 pairs done")

# ---------------------------------------------------------------------------
# Print
# ---------------------------------------------------------------------------
print()
hdr = (f"{'Pair':<10} {'AMBvac':>9} {'rAMB':>5}"
       f" {'AMB+SA':>9} {'rSA':>5}"
       f" {'CHARvac':>9} {'rCHAR':>5}"
       f" {'CHAR+SA':>9} {'rCSA':>5}"
       f" {'Upside':>9} {'rU':>5}")
print(hdr); print('-'*len(hdr))
for row in rows:
    def f(v): return f"{v:7.3f}" if not np.isnan(v) else "    nan"
    print(f"{row['Pair']:<10}"
          f" {row['AMBvac']:9.3f} {f(row['rAMB'])}"
          f" {row['AMB+SA']:9.3f} {f(row['rSA'])}"
          f" {row['CHARvac']:9.3f} {f(row['rCHAR'])}"
          f" {row['CHAR+SA']:9.3f} {f(row['rCSA'])}"
          f" {row['Upside']:9.3f} {f(row['rU'])}")

non_gly = [r for r in rows if 'GLY' not in r['Pair']]
for col in ['AMBvac','AMB+SA','CHARvac','CHAR+SA']:
    mad = np.mean([abs(r[col]-r['Upside']) for r in non_gly])
    bias = np.mean([r[col]-r['Upside'] for r in non_gly])
    print(f"  MAD({col:8s} vs Upside): {mad:.3f}  bias: {bias:+.3f} kcal/mol")

csv_path = "/tmp/sc_pair_comparison.csv"
flds = ['Pair','AMBvac','rAMB','AMB+SA','rSA','CHARvac','rCHAR','CHAR+SA','rCSA','Upside','rU']
with open(csv_path, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=flds)
    w.writeheader(); w.writerows(rows)
print(f"\nResults saved to {csv_path}")
