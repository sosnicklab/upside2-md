#!/usr/bin/env python
# Numpy oracle of the writers' index meaning (upside3.md 2.2).  For every node a config holds, the index datasets
# are rebuilt from first principles and compared, slot by slot as the C++ reader uses them, with what the writer
# stored.  The sources are the sequence; the atom layout of up.md 5 (N, CA, C at 3i, 3i+1, 3i+2; a hybrid's N, CA,
# C, O at 4i..4i+3, its BB beads at 4n+i, then the environment); the chain breaks the config records; the axis
# labels of the force-field tables the config names (/input/args, or for a hybrid the files martini_prepare_system
# fixes); and for a hybrid the dry-MARTINI topology.  Nothing from py/ is imported.  Conventions the C++ cannot
# define (the backbone classes of the coverage tables) are restated here and named as conventions.
# Each mismatch prints node, dataset, row, expected and stored (the first 200 rows of a dataset, then the count);
# a summary of every compared dataset follows, then the nodes the oracle does not check.
# usage: index_oracle.py CONFIG [CONFIG ...]
import sys, os, re
import numpy as np
import h5py
import tables as tb

UPSIDE_HOME = os.environ['UPSIDE_HOME']
PARAM = UPSIDE_HOME + '/parameters/'
ITP_DIR = UPSIDE_HOME + '/example/16.MARTINI/dryMARTINI_itp/'   # martini_prepare_system.py:1496
E_UP_KJ = 2.914952774272   # kJ/mol per E_up (CLAUDE.md unit table)
N_BIT_ROTAMER = 4          # bead_interaction.h:145
MAX_ROWS = 200

results = []   # (node, dataset, n_rows, n_bad)
notes = []

def show(x):
    return np.array2string(np.asarray(x), precision=6, separator=',', max_line_width=10000)

def compare(node, dset, expected, stored, atol=0.):
    e, s = np.asarray(expected), np.asarray(stored)
    if e.shape != s.shape:
        print('MISMATCH %s/%s shape: expected %s, stored %s' % (node, dset, e.shape, s.shape))
        results.append((node, dset, len(e) if e.ndim else 1, 'shape %s vs %s' % (e.shape, s.shape)))
        return False
    e2, s2 = e.reshape(len(e), -1) if e.ndim else e.reshape(1, 1), s.reshape(len(s), -1) if s.ndim else s.reshape(1, 1)
    if e2.dtype.kind in 'fc' or s2.dtype.kind in 'fc':
        bad = ~np.all(np.isclose(e2.astype('f8'), s2.astype('f8'), rtol=0., atol=atol, equal_nan=True), axis=1)
    else:
        bad = np.any(e2 != s2, axis=1)
    rows = np.nonzero(bad)[0]
    for r in rows[:MAX_ROWS]:
        print('MISMATCH %s/%s row %d: expected %s, stored %s' % (node, dset, r, show(e[r] if e.ndim else e),
                                                                  show(s[r] if s.ndim else s)))
    if len(rows) > MAX_ROWS:
        print('MISMATCH %s/%s: %d rows differ in all' % (node, dset, len(rows)))
    results.append((node, dset, len(e2), len(rows)))
    return not len(rows)

def note(text):
    print('NOTE ' + text)
    notes.append(text)

# ------------------------------------------------------------------------------------------------------------------
# config and layout
# ------------------------------------------------------------------------------------------------------------------

class Layout:
    # residue -> atom map of up.md 5, and the chains the config records
    def __init__(self, f):
        inp = f['input']
        self.seq_cpr = [x.decode().strip() for x in inp['sequence'][:]]
        self.seq = ['PRO' if x == 'CPR' else x for x in self.seq_cpr]
        self.n_res = n = len(self.seq)
        self.n_atom = inp['pos'].shape[0]
        self.hybrid = 'hybrid_control' in inp          # the engine's own switch (martini_hybrid.cpp:324)
        self.stride = 4 if self.hybrid else 3
        if 'chain_break' in inp:
            starts = inp['chain_break/chain_first_residue'][:]
        elif self.hybrid:
            cid = [x.decode() if isinstance(x, bytes) else str(x) for x in inp['hybrid_bb_map/bb_chain_id'][:]]
            starts = [i for i in range(1, n) if cid[i] != cid[i-1]]
        else:
            starts = []
        self.chain_starts = np.array(sorted(int(x) for x in starts), dtype=int)
        self.chain = np.searchsorted(self.chain_starts, np.arange(n), side='right')
        self.first = np.r_[True, self.chain[1:] != self.chain[:-1]]
        self.last = np.r_[self.chain[1:] != self.chain[:-1], True]
        r = np.arange(n)
        self.N, self.CA, self.C = self.stride*r, self.stride*r+1, self.stride*r+2
        if self.hybrid:
            self.O, self.BB = 4*r+3, 4*n+r

def read_args(path, group='args'):
    # the options upside_config (args) and advanced_config (adv_args) record
    with tb.open_file(path) as t:
        if group not in t.root.input: return {}
        a = t.root.input._f_get_child(group)._v_attrs
        return dict((k, a[k]) for k in a._v_attrnamesuser)

def arg_list(args, key):
    v = args.get(key)
    return set() if v is None else set(int(x) for x in np.atleast_1d(v))

def path_of(p):
    # option paths are recorded as given, relative to where the config was built (the repository root)
    p = str(p)
    return p if os.path.exists(p) else os.path.join(UPSIDE_HOME, p)

def libraries(args, hybrid):
    # the tables whose axes the index arrays point into
    if hybrid:   # martini_prepare_system.py:1498-1501, martini_prepare_system_lib.py:2247; martini.h5 per up.md 4
        return dict(sidechain=PARAM+'ff_2.1/sidechain.h5', rama=PARAM+'common/rama.dat', sheet=PARAM+'ff_2.1/sheet',
                    hbond=PARAM+'ff_2.1/hbond.h5', environment=None, membrane=None, martini=PARAM+'ff_2.1/martini.h5')
    s = lambda k: path_of(args[k]) if args.get(k) else None
    return dict(sidechain=s('rotamer_placement'), interaction=s('rotamer_interaction'), rama=s('rama_library'),
                sheet=s('rama_sheet_mixing_energy'), hbond=s('hbond_energy'), environment=s('environment_potential'),
                membrane=s('membrane_potential') or s('channel_membrane_potential'), martini=None)

def read_table(path):
    # the whitespace tables of the options: a header line, then rows
    rows = [l.split() for l in open(path_of(path)) if l.strip()]
    return [x.lower() for x in rows[0]], rows[1:]

def residue_list(s):
    # '2-10,12' (advanced_config's residue lists)
    out = []
    for part in s.split(','):
        a = part.split('-')
        out += list(range(int(a[0]), int(a[-1]) + 1))
    return out

def switching_prolines(args):
    # --trans-cis: the residues whose preceding peptide bond switches (the table's first column)
    return sorted(int(x[0]) for x in read_table(args['trans_cis'])[1]) if args.get('trans_cis') else []

def strings(ds):
    return [x.decode().strip() if isinstance(x, bytes) else str(x).strip() for x in ds[:]]

# ------------------------------------------------------------------------------------------------------------------
# backbone: bonded terms, dihedrals, frames
# ------------------------------------------------------------------------------------------------------------------

def check_backbone(pot, L, args):
    r = np.arange(L.n_res)
    inner = r[~L.last]
    if 'Distance3D' in pot:
        # coord_basic.cpp:33,42-43: id rows are atom pairs; Spring_bond reads row id of Distance3D (spring.cpp:39,55)
        rows, kind = [], []
        for i in r:
            rows += [(L.N[i], L.CA[i]), (L.CA[i], L.C[i])]; kind += [0, 1]
            if not L.last[i]: rows.append((L.C[i], L.N[i+1])); kind.append(2)
        rows, kind = np.array(rows), np.array(kind)
        stored = pot['Distance3D/id'][:]
        compare('Distance3D', 'id', np.sort(rows, axis=1), np.sort(stored, axis=1) if stored.shape == rows.shape else stored)
        if 'Spring_bond' in pot:
            sp = pot['Spring_bond']
            compare('Spring_bond', 'id', np.arange(len(stored)), sp['id'][:])
            if stored.shape == rows.shape: check_equil_by_kind('Spring_bond', sp['equil_dist'][:], kind, 'N-CA CA-C C-N')
            check_stiffness('Spring_bond', sp, args, 'bond_stiffness')
    if 'Angle' in pot:
        # coord_basic.cpp:222-225: the angle is at the THIRD atom (x1 = atom1-atom3, x2 = atom2-atom3)
        rows, kind = [], []
        for i in r:
            rows.append((L.N[i], L.C[i], L.CA[i])); kind.append(0)
            if not L.last[i]:
                rows += [(L.CA[i], L.N[i+1], L.C[i]), (L.C[i], L.CA[i+1], L.N[i+1])]; kind += [1, 2]
        rows, kind = np.array(rows), np.array(kind)
        stored = pot['Angle/id'][:]
        canon = lambda a: np.column_stack((np.sort(a[:, :2], axis=1), a[:, 2]))
        compare('Angle', 'id (ends sorted, vertex last)', canon(rows), canon(stored) if stored.shape == rows.shape else stored)
        if 'Spring_angle' in pot:
            sp = pot['Spring_angle']
            compare('Spring_angle', 'id', np.arange(len(stored)), sp['id'][:])
            if stored.shape == rows.shape: check_equil_by_kind('Spring_angle', sp['equil_dist'][:], kind, 'at-CA at-C at-N')
            check_stiffness('Spring_angle', sp, args, 'angle_stiffness')
    if 'Dihedral_phi' in pot:
        # coord_basic.cpp:327-330: atom -1 in slot 0 or 3 marks the missing angle of a chain end
        phi = np.column_stack((np.where(L.first, -1, L.C - L.stride), L.N, L.CA, L.C))
        compare('Dihedral_phi', 'id', phi, pot['Dihedral_phi/id'][:])
    if 'Dihedral_psi' in pot:
        psi = np.column_stack((L.N, L.CA, L.C, np.where(L.last, -1, L.N + L.stride)))
        compare('Dihedral_psi', 'id', psi, pot['Dihedral_psi/id'][:])
    omega = lambda i: np.column_stack((L.CA[i], L.C[i], L.N[i+1], L.CA[i+1]))
    pro = switching_prolines(args)
    fixed = np.array([i for i in inner if i+1 not in pro], dtype=int)   # a switching proline's omega moves out
    if 'Dihedral_omega' in pot:
        compare('Dihedral_omega', 'id', omega(fixed), pot['Dihedral_omega/id'][:])
        if 'Spring_omega' in pot:
            sp = pot['Spring_omega']
            compare('Spring_omega', 'id', np.arange(len(fixed)), sp['id'][:])
            cis = np.array([L.seq_cpr[i+1] == 'CPR' for i in fixed])
            compare('Spring_omega', 'equil_dist (cis-Pro 0, trans pi)', np.where(cis, 0., np.pi), sp['equil_dist'][:], 1e-6)
            check_stiffness('Spring_omega', sp, args, 'omega_stiffness')
    if 'Dihedral_OmegaTransCis' in pot:
        # the omega before each switching proline, in the table's (sorted) order; every lambda node below reads row k
        # as the k-th proline (coord_map.cpp:504, sigmoid_energy.cpp:72, rama_map_pot.cpp:157)
        compare('Dihedral_OmegaTransCis', 'id', omega(np.array(pro) - 1), pot['Dihedral_OmegaTransCis/id'][:])
        k = np.arange(len(pro))
        for node in [g for g in pot if g.startswith('SigmoidEnergy_') or g.startswith('SigmoidCoord_')]:
            compare(node, 'id (row of Dihedral_OmegaTransCis)', k, pot[node + '/id'][:])
        for node in ('Add_lambda_trans', 'Multiply_lambda_cis'):
            if node in pot:
                for d in ('id_pos1', 'id_pos2', 'id1_out', 'id2_out'):
                    compare(node, d, k, pot[node + '/' + d][:])
    if 'rama_coord' in pot:
        # coord_map.cpp:96-97: row i pairs phi row index_pos1[i] with psi row index_pos2[i]; every reader of rama_coord
        # takes its row i as residue i (hbond.cpp:525, placement.cpp:67, rama_map_pot.cpp:72)
        compare('rama_coord', 'index_pos1', r, pot['rama_coord/index_pos1'][:])
        compare('rama_coord', 'index_pos2', r, pot['rama_coord/index_pos2'][:])
    if 'affine_alignment' in pot:
        # eig.cpp:301-304, 326-345: atoms (N, CA, C) are fitted to ref_geom's slots in that order
        compare('affine_alignment', 'atoms', np.column_stack((L.N, L.CA, L.C)), pot['affine_alignment/atoms'][:])
        g = pot['affine_alignment/ref_geom'][:].astype('f8')
        d = np.column_stack((np.linalg.norm(g[:,0]-g[:,1], axis=1), np.linalg.norm(g[:,2]-g[:,1], axis=1),
                             np.linalg.norm(g[:,0]-g[:,2], axis=1), np.abs(g.sum(axis=1)).max(axis=1)))
        # slot 0 is N (N-CA 1.45 A), slot 2 is C (CA-C 1.53 A), and the three are centred (the fit centres the atoms)
        ok = (d[:,0] > 1.40) & (d[:,0] < 1.50) & (d[:,1] > 1.50) & (d[:,1] < 1.56) & (d[:,3] < 1e-5)
        compare('affine_alignment', 'ref_geom slots (N-CA, CA-C, N-C, centred)', np.ones(L.n_res, bool), ok)
    if 'backbone_pairs' in pot:
        check_backbone_pairs(pot, L)

def check_equil_by_kind(node, equil, kind, names):
    # each bond or angle type carries one equilibrium value, so a row shifted against its type shows up
    names = names.split()
    for k in range(kind.max()+1):
        v = equil[kind == k]
        compare(node, 'equil_dist of %s rows (one value per type)' % names[k], np.full(len(v), v[0]), v)
        print('INFO %s %s equil_dist %.6f' % (node, names[k], v[0]))

def check_stiffness(node, sp, args, key):
    k = sp['spring_const'][:]
    compare(node, 'spring_const', np.full(len(k), args[key] if key in args else k[0]), k, 1e-12)

def check_backbone_pairs(pot, L):
    # backbone_steric.cpp:66-74: id is the affine_alignment row; sequence_id (required) feeds the exclusion
    # acceptable_backbone_pair, |s1-s2| > 1 (:32-35); n_atom counts ref_pos slots (N, CA, C, CB) used (:120-128)
    bp = pot['backbone_pairs']
    r = np.arange(L.n_res)
    compare('backbone_pairs', 'id', r, bp['id'][:])
    gly = np.array([x == 'GLY' for x in L.seq])
    compare('backbone_pairs', 'n_atom', np.where(gly, 3, 4), bp['n_atom'][:])
    if compare('backbone_pairs', 'sequence_id present', True, 'sequence_id' in bp):
        s = bp['sequence_id'][:]
        expect_excl = (L.chain[:, None] == L.chain[None, :]) & (np.abs(r[:, None] - r[None, :]) <= 1)
        check_exclusion('backbone_pairs', 'sequence_id', expect_excl, np.abs(s[:, None] - s[None, :]) <= 1, r, r)
    ref = bp['ref_pos'][:].astype('f8')
    if 'affine_alignment' in pot:
        g = pot['affine_alignment/ref_geom'][:].astype('f8')
        compare('backbone_pairs', 'ref_pos[:, :3] (N, CA, C as affine ref_geom)', g, ref[:, :3], 1e-6)
    compare('backbone_pairs', 'ref_pos[:, 3] is NaN exactly for GLY', gly, np.isnan(ref[:, 3]).any(axis=1))

def check_exclusion(node, dset, expected, engine, res1, res2):
    # the engine excludes on the stored sequence ids; the physics excludes sequence neighbours within a chain only
    bad = np.argwhere(expected != engine)
    for a, b in bad[:MAX_ROWS]:
        print('MISMATCH %s/%s exclusion of element pair (%d, %d), residues (%d, %d): expected %s, engine %s'
              % (node, dset, a, b, res1[a], res2[b], 'excluded' if expected[a, b] else 'interacting',
                 'excluded' if engine[a, b] else 'interacting'))
    if len(bad) > MAX_ROWS:
        print('MISMATCH %s/%s exclusion: %d element pairs differ in all' % (node, dset, len(bad)))
    results.append((node, dset + ' exclusion', expected.size, len(bad)))

def check_pivots(f, pot, L):
    inp = f['input']
    r = np.arange(L.n_res)
    if 'pivot_moves' in inp:
        # monte_carlo_sampler.cpp:108-116: pivot_atom is (prev C, N, CA, C, next N); :147 rotates atoms in
        # [pivot_range[0], pivot_range[1]); :94-97 pivot_restype is the layer of proposal_pot
        pv = inp['pivot_moves']
        mid = r[~L.first & ~L.last]
        compare('pivot_moves', 'pivot_atom', np.column_stack((L.C[mid]-L.stride, L.N[mid], L.CA[mid], L.C[mid],
                                                             L.N[mid]+L.stride)), pv['pivot_atom'][:])
        compare('pivot_moves', 'pivot_range', np.column_stack((L.N[mid]+L.stride+1, np.full(len(mid), L.n_atom))),
                pv['pivot_range'][:])
        compare('pivot_moves', 'pivot_restype (map of the residue itself)', mid, pv['pivot_restype'][:])
        if 'rama_map_pot' in pot:
            compare('pivot_moves', 'proposal_pot (= rama_map_pot/rama_pot)', pot['rama_map_pot/rama_pot'][:],
                    pv['proposal_pot'][:], 1e-12)
    if 'jump_moves' in inp:
        # monte_carlo_sampler.cpp:188-194, 221: atom_range is [first atom, next chain's first atom)
        starts = np.r_[0, L.chain_starts, L.n_res] * L.stride
        cb = inp['chain_break']
        if 'rl_chains' in cb:
            k = int(cb['rl_chains'][0])
            expect = np.array([[0, starts[k]], [starts[k], starts[-1]]])
        else:
            expect = np.column_stack((starts[:-1], starts[1:]))
        compare('jump_moves', 'atom_range', expect, inp['jump_moves/atom_range'][:])

# ------------------------------------------------------------------------------------------------------------------
# Ramachandran maps
# ------------------------------------------------------------------------------------------------------------------

def mixture(weights, pots):
    # -log of the weight-normalised sum of exp(-V) (up.md 2.8a steps 1 and 3)
    w = np.asarray(weights, dtype='f8'); w = w / w.sum(axis=0)
    w = w.reshape(w.shape + (1,)*(np.ndim(pots) - w.ndim))
    return -np.log((w*np.exp(-np.asarray(pots, dtype='f8'))).sum(axis=0))

def library_side_maps(L, seq_cpr, group, allow_cpr, mode):
    names = [x.decode() if isinstance(x, bytes) else x for x in group._v_attrs.restype]
    dirs = [x.decode() if isinstance(x, bytes) else x for x in group._v_attrs.dir]
    idx = dict((x, i) for i, x in enumerate(names))
    pot, wt = group.dimer_pot[:].astype('f8'), group.dimer_weight[:].astype('f8')
    centre = lambda x: idx[x if (x != 'CPR' or allow_cpr) else 'PRO']
    neigh = lambda x: idx['PRO' if x == 'CPR' else x]   # a neighbour cis-proline reads as PRO (up.md 2.8)
    maps, weights = [], []
    for i, x in enumerate(seq_cpr):
        # left is the N-terminal neighbour i-1, right the C-terminal i+1; a chain end has one side only
        sides = ([] if L.first[i] else [('left', seq_cpr[i-1])]) + ([] if L.last[i] else [('right', seq_cpr[i+1])])
        V = [pot[centre(x), dirs.index(d), neigh(n)] for d, n in sides]
        W = [wt[centre(x), dirs.index(d), neigh(n)] for d, n in sides]
        if len(V) == 1: m = V[0]
        elif mode == 'product': m = V[0] + V[1] - pot[centre(x), dirs.index('right'), idx['ALL']]
        else: m = mixture(W, V)
        maps.append(m + np.log(np.exp(-m).sum()))
        weights.append(np.mean(W))
    return np.array(maps), np.array(weights), names

def secstr_basins(n_bin):
    # Convention of --secstr-bias (upside_config.py:788-802, A155): a window is 1 on an arc, its two 10 degree edge
    # sigmoids summed over the arc's periodic images; both basins take phi on (-240, 0) degrees, phi < 0 continued
    # past +-180 to +120, the helix psi on (-100, 50) and the sheet the rest of psi
    phi = np.linspace(-np.pi, np.pi, n_bin, endpoint=False)[:, None]
    psi = np.linspace(-np.pi, np.pi, n_bin, endpoint=False)[None, :]
    less = lambda a, b: 1./(1. + np.exp(-(b - a)/np.radians(10.)))
    arc = lambda x, lo, hi: sum(less(np.radians(lo), x + n*2.*np.pi) * less(x + n*2.*np.pi, np.radians(hi))
                                for n in (-1, 0, 1))
    helix_psi = arc(psi, -100., 50.)
    return dict(helix=arc(phi, -240., 0.) * helix_psi, sheet=arc(phi, -240., 0.) * (1. - helix_psi))

def rama_maps(L, lib, args, seq_cpr):
    # each residue's map as up.md 2.8a builds it: left/right mixture, normalisation, coil/sheet mixture, then the
    # Boltzmann-weighted shift; --secstr-bias adds its basin to the listed residues before the shift
    mode = str(args.get('rama_library_combining_rule') or 'mixture')
    seq = ['PRO' if x == 'CPR' else x for x in seq_cpr]
    with tb.open_file(lib['rama']) as t:
        coil, coil_w, _ = library_side_maps(L, seq_cpr, t.root.coil, True, mode)
        if lib['sheet']:
            # the sheet maps are always mixed (read_weighted_maps passes no mode; the sheet has no ALL neighbour)
            sheet, sheet_w, sheet_names = library_side_maps(L, seq_cpr, t.root.sheet, False, 'mixture')
            mix = np.loadtxt(lib['sheet'])
            maps = mixture([coil_w, sheet_w*np.exp(-mix[[sheet_names.index(x) for x in seq]])], [coil, sheet])
        else:
            maps = coil
    if args.get('secstr_bias'):
        basin = secstr_basins(maps.shape[1])
        for res, kind, energy in read_table(args['secstr_bias'])[1]:
            maps[int(res)] += float(energy) * basin[kind]
    return maps - (maps*np.exp(-maps)).sum(axis=(1, 2), keepdims=True)

def check_rama(pot, L, lib, args):
    r = np.arange(L.n_res)
    pro = switching_prolines(args)
    trans = [('PRO' if i in pro else x) for i, x in enumerate(L.seq_cpr)]
    plain = np.array([i for i in r if i not in pro], dtype=int)
    if 'rama_map_pot' in pot:
        # rama_map_pot.cpp:44-45, 72-76: residue_id is the rama_coord row, rama_map_id the layer of rama_pot, whose
        # axes are (phi, psi); a switching proline is left to RamaMap2_trans/cis
        g = pot['rama_map_pot']
        compare('rama_map_pot', 'residue_id', plain, g['residue_id'][:])
        compare('rama_map_pot', 'rama_map_id (one layer per residue)', plain, g['rama_map_id'][:])
        if 'rama_map_id_all' in g:
            compare('rama_map_pot', 'rama_map_id_all (read as pivot_restype)', r, g['rama_map_id_all'][:])
        if lib['rama']:
            compare('rama_map_pot', 'rama_pot[i] (residue i\'s left/right, coil/sheet map)', rama_maps(L, lib, args, trans),
                    g['rama_pot'][:], 2e-4)
    for node, state in (('RamaMap2_trans', 'PRO'), ('RamaMap2_cis', 'CPR')):
        if node in pot:
            # rama_map_pot.cpp:127-128, 148-157: row k is the k-th proline's rama_coord row and layer, weighted by row k
            # of the lambda node
            g = pot[node]
            compare(node, 'residue_id', np.array(pro), g['residue_id'][:])
            compare(node, 'rama_map_id', np.array(pro), g['rama_map_id'][:])
            seq = [(state if i in pro else x) for i, x in enumerate(trans)]
            compare(node, 'rama_pot[i] (the prolines %s)' % state, rama_maps(L, lib, args, seq), g['rama_pot'][:], 2e-4)
    if 'rama_map_pot_ref' in pot:
        g = pot['rama_map_pot_ref']
        compare('rama_map_pot_ref', 'residue_id', r, g['residue_id'][:])
        compare('rama_map_pot_ref', 'rama_map_id (one shared layer)', np.zeros(L.n_res, int), g['rama_map_id'][:])

# ------------------------------------------------------------------------------------------------------------------
# hydrogen bonds
# ------------------------------------------------------------------------------------------------------------------

def hbond_sites(L, args):
    # a donor H needs the previous residue's C and an NH (no proline); an acceptor O needs the next residue's N
    excl = arg_list(args, 'hbond_exclude_residues')
    r = np.arange(L.n_res)
    donors = np.array([i for i in r if not L.first[i] and L.seq[i] != 'PRO' and i not in excl], dtype=int)
    acceptors = np.array([i for i in r if not L.last[i] and i not in excl], dtype=int)
    return donors, acceptors

def check_hbond(pot, L, args, lib):
    donors, acceptors = hbond_sites(L, args)
    nd, na = len(donors), len(acceptors)
    if 'infer_H_O' in pot:
        # hbond.cpp:49-52, 77-86: H/O sits on atom[1]; atom[0] and atom[2] enter symmetrically; donors are rows
        # 0..nd-1 of the output, acceptors nd.. (:51)
        g = pot['infer_H_O']
        for side, res, rows in (('donors', donors, np.column_stack((L.C[donors]-L.stride, L.N[donors], L.CA[donors]))),
                                ('acceptors', acceptors, np.column_stack((L.CA[acceptors], L.C[acceptors],
                                                                          L.N[acceptors]+L.stride)))):
            st = g[side + '/id'][:]
            canon = lambda a: np.column_stack((a[:, 1], np.sort(a[:, [0, 2]], axis=1)))
            compare('infer_H_O', side + '/id (centre first, ends sorted)', canon(rows), canon(st) if st.shape == rows.shape else st)
            compare('infer_H_O', side + '/residue (read by the other writers)', res, g[side + '/residue'][:])
            bl = g[side + '/bond_length'][:]
            compare('infer_H_O', side + '/bond_length (one value per side)', np.full(len(bl), bl[0]), bl)
        print('INFO infer_H_O bond_length donors %.4f acceptors %.4f' % (g['donors/bond_length'][0], g['acceptors/bond_length'][0]))
        compare('infer_H_O', 'N-H shorter than C=O', True, g['donors/bond_length'][0] < g['acceptors/bond_length'][0])
    if 'protein_hbond' in pot:
        # hbond.cpp:329-340: output row nv copies infer_H_O row nv and an edge adds to row edge_indices1/2, so
        # index1/index2 must be the identity onto the donor and acceptor rows; id is the residue; sequence_id
        # excludes (:166-170)
        g = pot['protein_hbond']
        compare('protein_hbond', 'index1', np.arange(nd), g['index1'][:])
        compare('protein_hbond', 'index2', nd + np.arange(na), g['index2'][:])
        compare('protein_hbond', 'id1', donors, g['id1'][:])
        compare('protein_hbond', 'id2', acceptors, g['id2'][:])
        compare('protein_hbond', 'type1', np.zeros(nd, int), g['type1'][:])
        compare('protein_hbond', 'type2', np.zeros(na, int), g['type2'][:])
        igraph_exclusion('protein_hbond', g, L, donors, acceptors, 1)
    if 'hbond_energy' in pot:
        # hbond.cpp:469-472: rama_resid is read for its size only and rama_coord row i is residue i (:525); donor_resid
        # and acceptor_resid name each site's residue (:559-566), and a donor at r counts on peptide r-1, an acceptor
        # at r on peptide r (:571-577)
        g = pot['hbond_energy']
        compare('hbond_energy', 'donor_resid', donors, g['donor_resid'][:])
        compare('hbond_energy', 'acceptor_resid', acceptors, g['acceptor_resid'][:])
        compare('hbond_energy', 'rama_resid', np.arange(L.n_res), g['rama_resid'][:])
        if 'protein_hbond' in pot:
            compare('hbond_energy', 'len(donor_resid) = protein_hbond donors (only the sum is asserted, :484)',
                    len(pot['protein_hbond/index1']), len(g['donor_resid']))
        with tb.open_file(lib['hbond']) as t:
            classes = [x.decode() for x in t.root.class_restype[:]] if 'class_restype' in t.root else []
        if classes:
            compare('hbond_energy', 'residue_class', np.array([1+classes.index(x) if x in classes else 0 for x in L.seq]),
                    g['residue_class'][:])
        else:
            compare('hbond_energy', 'no residue_class without class_restype', False, 'residue_class' in g)

def igraph_exclusion(node, g, L, res1, res2, k):
    # interaction_graph.h:137 with the node's acceptable_id_pair on sequence_id1/2 (:490-492): elements whose
    # sequence ids differ by k or less never pair
    if not compare(node, 'sequence_id1 and sequence_id2 present', True, 'sequence_id1' in g and 'sequence_id2' in g):
        return
    s1, s2 = g['sequence_id1'][:], g['sequence_id2'][:]
    if len(s1) != len(res1) or len(s2) != len(res2): return
    expect = (L.chain[res1][:, None] == L.chain[res2][None, :]) & (np.abs(res1[:, None] - res2[None, :]) <= k)
    check_exclusion(node, 'sequence_id1/sequence_id2', expect, np.abs(s1[:, None] - s2[None, :]) <= k, res1, res2)

# ------------------------------------------------------------------------------------------------------------------
# side chains: placement, rotamer, coverage
# ------------------------------------------------------------------------------------------------------------------

def rotamer_rows(L, sclib, excluded, fixed={}):
    # one row per bead of every rotamer of every residue, in library order (rotamer-major); the rotamer node's id
    # packs (rotamer, n_rot, node), node counting the residues of each n_rot (rotamer.cpp:884-888).  A residue in
    # fixed keeps the one rotamer given.
    names = [x.decode() for x in sclib.restype_order[:]]
    start_stop = sclib.rotamer_start_stop_bead[:]
    rows, count = [], {}
    for i, x in enumerate(L.seq):
        if i in excluded: continue
        start, stop, nb = start_stop[names.index(x)]
        if i in fixed: start, stop = start + nb*fixed[i], start + nb*(fixed[i] + 1)
        n_rot = (stop - start) // nb
        node = count.get(n_rot, 0); count[n_rot] = node + 1
        for rel in range(stop - start):
            rows.append((i, start + rel, rel // nb, rel % nb, n_rot, node))
    return np.array(rows, dtype=int).reshape(-1, 6)

def radial_scaled(table, args):
    # --rot-scale multiplies the radial knots, which follow the two beads' 15 angular knots each (up.md 2.1)
    t = np.array(table, dtype='f8')
    t[..., 30:] *= float(args.get('rot_scale', 1.))
    return t

def fixed_rotamers(node, L, sclib, args, stored_res, stored_layer):
    # --fix-rotamer: the state each listed residue keeps is the writer's choice from its chi angles (not rebuilt
    # here); what is checked is that only listed residues are reduced, each to exactly one whole library rotamer
    names = [x.decode() for x in sclib.restype_order[:]]
    listed = set(int(x[0]) for x in read_table(args['fix_rotamer'])[1])
    fixed = {}
    for i in sorted(set(stored_res.tolist())):
        start, stop, nb = sclib.rotamer_start_stop_bead[:][names.index(L.seq[i])]
        layers = stored_layer[stored_res == i]
        if len(layers) == stop - start: continue
        ok = i in listed and len(layers) == nb and (layers[0] - start) % nb == 0 and np.all(np.diff(layers) == 1)
        compare(node, '--fix-rotamer: residue %d reduced to one whole rotamer of a listed residue' % i, True, bool(ok))
        fixed[i] = int((layers[0] - start) // nb)
    return fixed

def backbone_class(L):
    # Convention of the coverage tables (upside_config.py:52-56 backbone_group4, :1198-1213, :1240-1247): 0 generic,
    # 1 GLY, 2 PRO, 3 a residue bonded to a following proline (pre-proline takes precedence)
    cls = np.array([1 if x == 'GLY' else 2 if x == 'PRO' else 0 for x in L.seq])
    pre = np.array([not L.last[i] and L.seq[i+1] == 'PRO' for i in range(L.n_res)])
    return np.where(pre, 3, cls)

def check_sidechain(pot, L, args, lib):
    sc_node = 'placement_fixed_point_vector_only'
    if sc_node not in pot or not lib['sidechain']: return None
    sclib = tb.open_file(lib['sidechain'])
    excluded = arg_list(args, 'rotamer_exclude_residues')
    g = pot[sc_node]
    stored_res = g['affine_residue'][:]
    fixed = fixed_rotamers(sc_node, L, sclib.root, args, stored_res, g['layer_index'][:]) if args.get('fix_rotamer') else {}
    rows = rotamer_rows(L, sclib.root, excluded, fixed)
    compare(sc_node, 'rows per residue', np.bincount(rows[:, 0], minlength=L.n_res),
            np.bincount(stored_res, minlength=L.n_res))
    absent = sorted(set(range(L.n_res)) - set(stored_res.tolist()) - excluded)
    if absent:
        # check the rows that are there as what they claim to be
        note('%s has no rows for residues %s (%s); the other checks use the rows present'
             % (sc_node, absent, [L.seq[i] for i in absent]))
        rows = rotamer_rows(L, sclib.root, excluded | set(absent), fixed)
    res, layer, rot, bead, n_rot, node_id = rows.T
    bead_order = [x.decode() for x in sclib.root.bead_order[:]]
    btype = np.array([bead_order.index('%s_%d' % (L.seq[i], b)) for i, b in zip(res, bead)])

    # placement.cpp:252, 272, 129, 140: affine_residue is the affine_alignment row, layer_index the placement_data row
    compare(sc_node, 'affine_residue', res, stored_res)
    compare(sc_node, 'layer_index', layer, g['layer_index'][:])
    compare(sc_node, 'placement_data (= library rotamer_center_fixed, so layers are library rows)',
            sclib.root.rotamer_center_fixed[:, :6], g['placement_data'][:], 2e-6)
    pl_node = next((x for x in ('placement_scalar', 'placement_fixed_scalar') if x in pot), None)
    if pl_node:
        # rows must align with the bead rows: rotamer.cpp:715 sizes them, :889-892 reads both at loc1
        p = pot[pl_node]
        compare(pl_node, 'affine_residue', res, p['affine_residue'][:])
        compare(pl_node, 'layer_index', layer, p['layer_index'][:])
        if pl_node == 'placement_scalar':
            # placement.cpp:67-75: rama_residue is the rama_coord row; the table axes are (layer, phi, psi)
            compare(pl_node, 'rama_residue', res, p['rama_residue'][:])
            e = -np.log(sclib.root.rotamer_prob[:].transpose((2, 0, 1)))[..., None]
            compare(pl_node, 'placement_data (= -log library rotamer_prob, axes layer, phi, psi)', e,
                    p['placement_data'][:], 1e-5)
    if 'weighted_pos' in pot:
        # environment.cpp:350-353: position of row index_pos, weight exp(-1-body) of row index_weight; both the same
        # rotamer row
        w = pot['weighted_pos']
        compare('weighted_pos', 'index_pos', np.arange(len(rows)), w['index_pos'][:])
        compare('weighted_pos', 'index_weight', np.arange(len(rows)), w['index_weight'][:])
        compare('weighted_pos', 'arguments', [sc_node, pl_node], [x.decode() for x in w.attrs['arguments']])
    if 'rotamer' in pot:
        # interaction_graph.h:378-381: index is the placement row (also the row of every probability node,
        # rotamer.cpp:889), type the pair_interaction row; id decoded at rotamer.cpp:886-888; sequence_id excludes
        # a residue's own beads (bead_interaction.h:148-150)
        pi = pot['rotamer/pair_interaction']
        compare('rotamer', 'pair_interaction/index', np.arange(len(rows)), pi['index'][:])
        compare('rotamer', 'pair_interaction/type (bead_order)', btype, pi['type'][:])
        packed = rot + (n_rot << N_BIT_ROTAMER) + (node_id << 2*N_BIT_ROTAMER)
        compare('rotamer', 'pair_interaction/id (rotamer | n_rot<<4 | node<<8)', packed, pi['id'][:])
        compare('rotamer', 'n_rot in {1, 3, 6} (rotamer.cpp:632)', np.ones(len(rows), bool), np.isin(n_rot, (1, 3, 6)))
        with tb.open_file(lib.get('interaction') or lib['sidechain']) as t:
            compare('rotamer', 'pair_interaction/interaction_param (= library pair_interaction, --rot-scale applied)',
                    radial_scaled(t.root.pair_interaction[:], args), pi['interaction_param'][:], 2e-6)
        args_ = [x.decode() for x in pot['rotamer'].attrs['arguments']]
        compare('rotamer', 'arguments[0] (the bead positions)', sc_node, args_[0])
        if compare('rotamer', 'pair_interaction/sequence_id present', True, 'sequence_id' in pi):
            sid = pi['sequence_id'][:]
            if len(sid) == len(rows):
                expect = res[:, None] == res[None, :]
                check_exclusion('rotamer', 'pair_interaction/sequence_id', expect, sid[:, None] == sid[None, :], res, res)

    donors, acceptors = hbond_sites(L, args)
    nd, na = len(donors), len(acceptors)
    cls = backbone_class(L)
    if 'hbond_coverage' in pot:
        # hbond.cpp:383-399: group 1 the protein_hbond rows (types: coverage_interaction rows), group 2 the bead rows;
        # output row j is element j of group 2 and the rotamer reads it at the bead row, so index2 is the identity
        g = pot['hbond_coverage']
        sites = np.r_[donors, acceptors]
        compare('hbond_coverage', 'index1', np.arange(nd+na), g['index1'][:])
        compare('hbond_coverage', 'id1', sites, g['id1'][:])
        compare('hbond_coverage', 'type1 (2*backbone class + acceptor)', np.r_[2*cls[donors], 2*cls[acceptors]+1], g['type1'][:])
        compare('hbond_coverage', 'index2', np.arange(len(rows)), g['index2'][:])
        compare('hbond_coverage', 'id2', res, g['id2'][:])
        compare('hbond_coverage', 'type2 (bead_order)', btype, g['type2'][:])
        compare('hbond_coverage', 'interaction_param (= library coverage_interaction, --rot-scale applied)',
                radial_scaled(sclib.root.coverage_interaction[:], args), g['interaction_param'][:], 2e-6)
        igraph_exclusion('hbond_coverage', g, L, sites, res, 2)
    if 'placement_fixed_point_vector_scalar' in pot:
        g = pot['placement_fixed_point_vector_scalar']
        compare('placement_fixed_point_vector_scalar', 'affine_residue', np.repeat(np.arange(L.n_res), 3), g['affine_residue'][:])
        compare('placement_fixed_point_vector_scalar', 'layer_index', np.tile(np.arange(3), L.n_res), g['layer_index'][:])
        compare('placement_fixed_point_vector_scalar', 'placement_data (= library hydrophobe_placement)',
                sclib.root.hydrophobe_placement[:], g['placement_data'][:], 2e-6)
    if 'hbond_coverage_hydrophobe' in pot:
        g = pot['hbond_coverage_hydrophobe']
        bb_res = np.repeat(np.arange(L.n_res), 3)
        compare('hbond_coverage_hydrophobe', 'index1', np.arange(3*L.n_res), g['index1'][:])
        compare('hbond_coverage_hydrophobe', 'id1', bb_res, g['id1'][:])
        compare('hbond_coverage_hydrophobe', 'type1 (3*backbone class + hydrophobe layer)',
                3*cls[bb_res] + np.tile(np.arange(3), L.n_res), g['type1'][:])
        compare('hbond_coverage_hydrophobe', 'index2', np.arange(len(rows)), g['index2'][:])
        compare('hbond_coverage_hydrophobe', 'id2', res, g['id2'][:])
        compare('hbond_coverage_hydrophobe', 'type2 (bead_order)', btype, g['type2'][:])
        compare('hbond_coverage_hydrophobe', 'interaction_param (= library hydrophobe_interaction, --rot-scale applied)',
                radial_scaled(sclib.root.hydrophobe_interaction[:], args), g['interaction_param'][:], 2e-6)
        igraph_exclusion('hbond_coverage_hydrophobe', g, L, bb_res, res, 2)
    sclib.close()
    return dict(res=res, layer=layer, rot=rot, bead=bead, n_rot=n_rot, btype=btype, node=sc_node)

def cb_frame(pot):
    # CA and C of the residue frame: affine_alignment's ref_geom slots 1 and 2
    g = pot['affine_alignment/ref_geom'][:].astype('f8')
    return g[:, 1], g[:, 2]

def check_cb(pot, L, args):
    r = np.arange(L.n_res)
    for node in ('placement_fixed_point_vector_only_CB', 'placement_fixed_point_only_CB'):
        if node not in pot: continue
        g = pot[node]
        compare(node, 'affine_residue', r, g['affine_residue'][:])
        compare(node, 'layer_index', np.zeros(L.n_res, int), g['layer_index'][:])
        d = g['placement_data'][:].astype('f8')
        if 'backbone_pairs' in pot:
            cb = pot['backbone_pairs/ref_pos'][:, 3].astype('f8')
            cb = cb[np.isfinite(cb).all(axis=1)][0]
            compare(node, 'placement_data[:, :3] (= backbone_pairs CB)', cb[None], d[:, :3], 1e-6)
            if d.shape[1] == 6 and 'affine_alignment' in pot:
                ca, c = cb_frame(pot)
                origin = c[0] if args.get('vector_CA_CO') else ca[0]
                u = (cb - origin) / np.linalg.norm(cb - origin)
                compare(node, 'placement_data[:, 3:] (unit CA->CB, or C->CB with --vector-CA-CO)', u[None], d[:, 3:], 1e-6)

def check_environment(pot, L, args, lib, sc):
    if not lib['environment'] or sc is None: return
    with tb.open_file(lib['environment']) as t:
        names = [x.decode() for x in t.root.restype_order[:]]
        coverage_param = t.root.coverage_param[:]
    etype = np.array([names.index(x) for x in L.seq])
    excl = arg_list(args, 'env_exclude_residues')
    cb_res = np.array([i for i in range(L.n_res) if i not in excl])
    n_sc = len(sc['res'])
    if 'environment_coverage_sc' in pot:
        # environment.cpp:193-226: group 1 the CB rows (type: coverage_param row of the CB's residue type), group 2
        # the weighted_pos rows; output element e*20 + aa_types[id2]; aa_types is indexed by the covering residue
        g = pot['environment_coverage_sc']
        compare('environment_coverage_sc', 'index1 (CB row = residue)', cb_res, g['index1'][:])
        compare('environment_coverage_sc', 'id1', cb_res, g['id1'][:])
        compare('environment_coverage_sc', 'type1 (environment.h5 restype)', etype[cb_res], g['type1'][:])
        compare('environment_coverage_sc', 'index2', np.arange(n_sc), g['index2'][:])
        compare('environment_coverage_sc', 'id2', sc['res'], g['id2'][:])
        compare('environment_coverage_sc', 'type2', np.zeros(n_sc, int), g['type2'][:])
        compare('environment_coverage_sc', 'aa_types (residue type by residue)', etype, g['aa_types'][:])
        compare('environment_coverage_sc', 'interaction_param (= environment.h5 coverage_param)', coverage_param,
                g['interaction_param'][:], 1e-12)
        igraph_exclusion('environment_coverage_sc', g, L, cb_res, sc['res'], 2)
    for node in ('sigmoid_coupling_environment', 'nonlinear_coupling_environment'):
        if node in pot:
            # environment.cpp:634-641 (sigmoid): element nr of the coverage takes coupling_types[nr]'s parameters
            compare(node, 'coupling_types (CB residue type by element)', etype[cb_res], pot[node + '/coupling_types'][:])
    donors, acceptors = hbond_sites(L, args)
    nd, na = len(donors), len(acceptors)
    sites = np.r_[donors, acceptors]
    if 'environment_coverage_hb' in pot:
        g = pot['environment_coverage_hb']
        compare('environment_coverage_hb', 'index1 (infer_H_O rows)', np.arange(nd+na), g['index1'][:])
        compare('environment_coverage_hb', 'id1', sites, g['id1'][:])
        compare('environment_coverage_hb', 'type1', np.zeros(nd+na, int), g['type1'][:])
        compare('environment_coverage_hb', 'index2', np.arange(n_sc), g['index2'][:])
        compare('environment_coverage_hb', 'id2', sc['res'], g['id2'][:])
        compare('environment_coverage_hb', 'type2', np.zeros(n_sc, int), g['type2'][:])
        compare('environment_coverage_hb', 'aa_types (residue type by residue)', etype, g['aa_types'][:])
        igraph_exclusion('environment_coverage_hb', g, L, sites, sc['res'], 2)
    for node, i1, r1, i2, r2 in (('hb_environment_coverage_hn', np.arange(nd), donors, nd+np.arange(na), acceptors),
                                 ('hb_environment_coverage_oc', nd+np.arange(na), acceptors, np.arange(nd), donors)):
        # environment.cpp:252-267: both groups are infer_H_O rows; output element e is group-1 element e
        if node not in pot: continue
        g = pot[node]
        compare(node, 'index1', i1, g['index1'][:]); compare(node, 'id1', r1, g['id1'][:])
        compare(node, 'index2', i2, g['index2'][:]); compare(node, 'id2', r2, g['id2'][:])
        compare(node, 'type1', np.zeros(len(i1), int), g['type1'][:]); compare(node, 'type2', np.zeros(len(i2), int), g['type2'][:])
        igraph_exclusion(node, g, L, r1, r2, 2)
    if 'cat_pos_bb_coverage' in pot:
        # coord_map.cpp:130-131, 137-144: rows of hn then rows of oc
        g = pot['cat_pos_bb_coverage']
        compare('cat_pos_bb_coverage', 'index_pos1', np.arange(nd), g['index_pos1'][:])
        compare('cat_pos_bb_coverage', 'index_pos2', np.arange(na), g['index_pos2'][:])
    if all(x in pot for x in ('bb_sigmoid_coupling_environment', 'environment_coverage_hb', 'cat_pos_bb_coverage',
                              'hb_environment_coverage_hn', 'hb_environment_coverage_oc')):
        # environment.cpp:739-768: the first argument is the side-chain coverage (element nr, 20 types), the second the
        # H-bond coverage (element nr); both elements nr must be the same backbone site
        a = [x.decode() for x in pot['bb_sigmoid_coupling_environment'].attrs['arguments']]
        compare('bb_sigmoid_coupling_environment', 'arguments', ['environment_coverage_hb', 'cat_pos_bb_coverage'], a)
        ip1, ip2 = pot['cat_pos_bb_coverage/index_pos1'][:], pot['cat_pos_bb_coverage/index_pos2'][:]
        via_cat = np.r_[pot['hb_environment_coverage_hn/index1'][:][ip1], pot['hb_environment_coverage_oc/index1'][:][ip2]]
        compare('bb_sigmoid_coupling_environment', 'site of element nr: environment_coverage_hb vs cat_pos_bb_coverage',
                pot['environment_coverage_hb/index1'][:], via_cat)

def check_membrane(pot, L, args, lib):
    donors, acceptors = hbond_sites(L, args)
    sites = np.r_[donors, acceptors]
    excl = arg_list(args, 'membrane_exclude_residues')
    res = np.array([i for i in range(L.n_res) if i not in excl])
    names = None
    if lib['membrane']:
        with tb.open_file(lib['membrane']) as t:
            names = [x.decode() for x in t.root.names[:]]
    for node, tdset in (('cb_membrane_potential', 'res_type'), ('cb_surf_membrane_potential', 'res_type')):
        if node not in pot or names is None: continue
        # membrane_potential.cpp:283-284, 316-322, 343 (surface variant :744-745, 783-788, 840): cb_index is the CB
        # row of the first argument, the element of environment_coverage_sc
        # (env_cov(0, ri*20+aa)) and the surface row; the residue type is the membrane.h5 row
        g = pot[node]
        cb_index = g['cb_index'][:].astype(int)
        compare(node, 'cb_index', res, cb_index)
        compare(node, tdset + ' (membrane.h5 names)', np.array([names.index(L.seq[i]) for i in res]), g[tdset][:])
        a = [x.decode() for x in g.attrs['arguments']]
        if a[1] == 'environment_coverage_sc' and 'environment_coverage_sc' in pot and len(cb_index):
            id1 = pot['environment_coverage_sc/id1'][:]
            cb_res = pot[a[0] + '/affine_residue'][:]
            ok = cb_index < len(id1)
            compare(node, 'residue of cb_index in the CB node vs in environment_coverage_sc', cb_res[cb_index[ok]],
                    id1[cb_index[ok]])
    rows = np.array([i for i, x in enumerate(sites) if x not in excl], dtype=int)
    for node in ('hb_membrane_potential', 'hb_surf_membrane_potential'):
        if node not in pot: continue
        # membrane_potential.cpp:529-530, 549-552 (surface variant :1050-1052, 1083-1105): hb_index is the
        # protein_hbond row, hb_type[nv] (0 donor, 1 acceptor) selects the coeff pair of that slot, hb_res the residue
        # (coverage element and surface row); an excluded residue's H and O leave with its CB (A168)
        g = pot[node]
        compare(node, 'hb_index', rows, g['hb_index'][:].astype(int))
        compare(node, 'hb_type (by slot)', (rows >= len(donors)).astype(int), g['hb_type'][:])
        if 'hb_res' in g:
            compare(node, 'hb_res', sites[rows], g['hb_res'][:])
    if 'cat_pos_bb_cb' in pot:
        # surface.cpp:185-219: row 3nr+j of its input is residue nr's N, CA, C and row index_cb[nr] its CB
        g = pot['cat_pos_bb_cb']
        compare('cat_pos_bb_cb', 'index_pos1 (N, CA, C of each residue)', np.column_stack((L.N, L.CA, L.C)).ravel(),
                g['index_pos1'][:])
        compare('cat_pos_bb_cb', 'index_pos2 (CB rows)', np.arange(L.n_res), g['index_pos2'][:])
    if 'surface' in pot:
        g = pot['surface']
        compare('surface', 'index_cb (cat_pos_bb_cb row of each CB)', 3*L.n_res + np.arange(L.n_res), g['index_cb'][:])
        inc = arg_list(args, 'surface_included_residues')
        compare('surface', 'exclusion (1 with an included list)', int(bool(inc)), int(g.attrs['exclusion']))
        if inc:
            compare('surface', 'included_list', np.array(sorted(inc)), g['included_list'][:])
    if 'hbbb_coverage' in pot:
        # environment.cpp:286-305 with HbondEnvironmentCoverageInteraction2: group 1 the donor rows of infer_H_O,
        # group 2 the backbone atoms of pos; ids are residues (exclusion |sequence_id1-sequence_id2| > 2)
        g = pot['hbbb_coverage']
        atoms = np.column_stack((L.N, L.CA, L.C)).ravel()
        compare('hbbb_coverage', 'index1', np.arange(len(donors)), g['index1'][:])
        compare('hbbb_coverage', 'id1', donors, g['id1'][:])
        compare('hbbb_coverage', 'index2 (backbone atoms)', atoms, g['index2'][:])
        compare('hbbb_coverage', 'id2 (residue of the atom)', np.repeat(np.arange(L.n_res), 3), g['id2'][:])
        igraph_exclusion('hbbb_coverage', g, L, donors, np.repeat(np.arange(L.n_res), 3), 2)

# ------------------------------------------------------------------------------------------------------------------
# advanced_config: restraints, groups, pulling, contacts
# ------------------------------------------------------------------------------------------------------------------

AXES = dict(x=[0], y=[1], z=[2], xy=[0, 1], xz=[0, 2], yz=[1, 2], xyz=[0, 1, 2])

def distance_node(pot, name, dims):
    # Distance1D/2D/3D_<name> by the table's axes; coord_basic.cpp:33 (3D), :92 (2D), :150 (1D) read id rows of
    # (row of argument 1, row of argument 2) and dim1/dim2 select the axes (:83-84, :142)
    node = 'Distance%dD_%s' % (len(dims), name)
    if node not in pot:
        compare(name, 'distance node %s present' % node, True, False)
        return None
    g = pot[node]
    attrs = [int(g.attrs['dim%d' % (k+1)]) for k in range(len(dims))] if len(dims) < 3 else dims
    compare(node, 'dim1/dim2 (the table\'s axes)', np.array(dims), np.array(attrs))
    return g

def check_spring(pot, name, kind, n, equil, k, wall_type=None):
    # spring.cpp:39-41, 55 (Spring) and :111-114, 126-131 (WallSpring): row nt acts on row id[nt] of the distance node
    node = '%s_%s' % (kind, name)
    g = pot[node]
    compare(node, 'id', np.arange(n), g['id'][:])
    compare(node, 'equil_dist', np.array(equil, dtype='f8'), g['equil_dist'][:], 1e-5)
    compare(node, 'spring_const', np.array(k, dtype='f8'), g['spring_const'][:], 1e-9)
    if wall_type is not None:
        compare(node, 'wall_type (0 inner wall, 1 outer wall)', np.array(wall_type, dtype='f8'), g['wall_type'][:])

def const_rows(pot, node, rows, points, dims):
    # a fixed point is the Const3D row named in the distance's second slot (coord_const.cpp:41-43); on the axes the
    # table gives it carries the table's value, elsewhere 0
    value = pot[node + '/value'][:]
    want = np.zeros((len(points), 3)); want[:, dims] = np.array(points, dtype='f8').reshape(len(points), -1)
    compare(node, 'value[id] (the table\'s reference point)', want, value[np.asarray(rows, dtype=int)], 1e-9)

def check_advanced(f, pot, L, adv):
    CA = L.CA
    for node in [g for g in pot if g.startswith('Const3D')]:
        # coord_const.cpp:13-16: the node's size is id's, its rows are value's, so the two must agree
        compare(node, 'id (one per value row)', np.arange(len(pot[node + '/value'])), pot[node + '/id'][:])
    for opt, name, kind, n_fix in (('fixed_wall', 'wall', 'WallSpring', 4), ('pair_wall', 'pair_wall', 'WallSpring', 5),
                                   ('fixed_spring', 'spring', 'Spring', 3), ('pair_spring', 'pair_spring', 'Spring', 4),
                                   ('nail_spring', 'nail', 'Spring', 2)):
        if not adv.get(opt): continue
        head, table = read_table(adv[opt])
        pair = opt.startswith('pair')
        rows = table if pair else [r for r in table for _ in residue_list(r[0])]
        res = None if pair else sum([residue_list(r[0]) for r in table], [])
        dims = [['x0', 'y0', 'z0'].index(x) for x in head[n_fix:]] if opt in ('fixed_wall', 'fixed_spring') else AXES[head[-1]]
        g = distance_node(pot, name, dims)
        if g is None: continue
        node, ids = g.name.split('/')[-1], g['id'][:].astype(int)
        if pair:
            compare(node, 'id (CA of residue1, CA of residue2)', CA[[[int(r[0]), int(r[1])] for r in rows]], ids)
        else:
            compare(node, 'id[:, 0] (CA of each residue)', CA[res], ids[:, 0])
            if opt == 'nail_spring':
                # Const3D_pos row k holds the config's position of the CA it nails (coord_const.cpp:13-16): the nail
                # ties the CA to its own starting position
                compare('Const3D_pos', 'value[slot 1 of %s] (the same CA at the start)' % node,
                        f['input/pos'][:, :, 0][CA[res]], pot['Const3D_pos/value'][:][ids[:, 1]], 1e-6)
            else:
                const_rows(pot, 'Const3D', ids[:, 1], [[float(v) for v in r[n_fix:]] for r in rows], dims)
        c = 2 if pair else 1   # the radius column
        if opt == 'nail_spring':
            check_spring(pot, name, kind, len(rows), [0.]*len(rows), [float(r[1]) for r in rows])
        else:
            check_spring(pot, name, kind, len(rows), [float(r[c]) for r in rows], [float(r[c+1]) for r in rows],
                         [float(r[c+2]) for r in rows] if kind == 'WallSpring' else None)

    if adv.get('group'):
        # coord_basic.cpp:399-401, 423-424: group g averages pos rows index_pos[border[g]:border[g+1]]
        groups = [residue_list(r[0]) for r in read_table(adv['group'])[1]]
        g = pot['GroupCenter_xyz']
        compare('GroupCenter_xyz', 'index_pos (CA of each group\'s residues)', CA[sum(groups, [])], g['index_pos'][:])
        compare('GroupCenter_xyz', 'border', np.cumsum([0] + [len(x) for x in groups]), g['border'][:])
        compare('GroupCenter_xyz', 'index_dim', np.arange(3), g['index_dim'][:])
    for opt, name, n_fix, kind in (('group_dist_spring', 'group_dist', 2, 'pair'), ('group_nail_spring', 'group_nail', 1, 'nail'),
                                   ('group_const_spring', 'group_const', 1, 'const')):
        if not adv.get(opt): continue
        head, rows = read_table(adv[opt])
        dims = [['x0', 'y0', 'z0'].index(x) for x in head[3:]] if kind == 'const' else AXES[head[-1]]
        g = distance_node(pot, name, dims)
        if g is None: continue
        ids = g['id'][:]
        nm = g.name.split('/')[-1]
        if kind == 'pair':
            compare(nm, 'id (rows of GroupCenter_xyz)', np.array([[int(r[0]), int(r[1])] for r in rows]), ids.astype(int))
            check_spring(pot, name, 'Spring', len(rows), [float(r[2]) for r in rows], [float(r[3]) for r in rows])
        elif kind == 'nail':
            compare(nm, 'id[:, 0] (row of GroupCenter_xyz)', np.array([int(r[0]) for r in rows]), ids[:, 0].astype(int))
            compare('Const3D_GroupCenter_xyz', 'id[slot 1 of %s] (the same group)' % nm, np.array([int(r[0]) for r in rows]),
                    pot['Const3D_GroupCenter_xyz/id'][:][ids[:, 1].astype(int)])
            check_spring(pot, name, 'Spring', len(rows), [0.]*len(rows), [float(r[1]) for r in rows])
        else:
            compare(nm, 'id[:, 0] (row of GroupCenter_xyz)', np.array([int(r[0]) for r in rows]), ids[:, 0].astype(int))
            const_rows(pot, 'Const3D', ids[:, 1], [[float(v) for v in r[3:]] for r in rows], dims)
            check_spring(pot, name, 'Spring', len(rows), [float(r[1]) for r in rows], [float(r[2]) for r in rows])
        if ids.dtype.kind == 'f':
            note('%s/id is stored as float (read through an int conversion)' % nm)

    groups = [residue_list(x) for x in np.atleast_1d(adv.get('restraint_group', []))]
    if adv.get('apply_restraint_group_to_each_chain'):
        # one group per chain, or the receptor chains and the ligand chains as two groups
        starts = np.r_[0, L.chain_starts, L.n_res]
        cb = f['input/chain_break']
        k = int(cb['rl_chains'][0]) if 'rl_chains' in cb else None
        cuts = [(0, k), (k, len(starts)-1)] if k else [(c, c+1) for c in range(len(starts)-1)]
        groups += [list(range(starts[a], starts[b])) for a, b in cuts]
    x = f['input/pos'][:, :, 0].astype('f8')
    for k, grp in enumerate(groups):
        # make_restraint_group: springs between backbone atoms of different residues of the group, at their
        # distance in /input/pos
        node = 'Distance3D_restraint%d' % k
        ids = pot[node + '/id'][:].astype(int)
        atoms = np.column_stack((L.N, L.CA, L.C))[grp].ravel()
        res_of = np.full(L.n_atom, -1); res_of[np.column_stack((L.N, L.CA, L.C)).ravel()] = np.repeat(np.arange(L.n_res), 3)
        ok = np.isin(ids, atoms).all(axis=1) & (ids[:, 0] < ids[:, 1]) & (res_of[ids[:, 0]] != res_of[ids[:, 1]])
        compare(node, 'id: backbone atoms of the group, ordered, of different residues', np.ones(len(ids), bool), ok)
        compare(node, 'id rows distinct', len(ids), len(set(map(tuple, ids.tolist()))))
        check_spring(pot, 'restraint%d' % k, 'Spring', len(ids), np.linalg.norm(x[ids[:, 0]] - x[ids[:, 1]], axis=1),
                     np.full(len(ids), float(adv['restraint_spring_constant'])))

    if 'Distance3D_cavity' in pot:
        # write_cavity_radial: every atom against the Const3D origin, an outer wall at the cavity radius
        ids = pot['Distance3D_cavity/id'][:].astype(int)
        compare('Distance3D_cavity', 'id[:, 0] (every atom)', np.arange(L.n_atom), ids[:, 0])
        const_rows(pot, 'Const3D', ids[:, 1], [[0., 0., 0.]]*len(ids), [0, 1, 2])
        radius = pot['WallSpring_cavity_radial/equil_dist'][:]
        if adv.get('cavity_radius_from_config') and os.path.exists(path_of(adv['cavity_radius_from_config'])):
            with h5py.File(path_of(adv['cavity_radius_from_config']), 'r') as other:
                radius = np.full(len(ids), other['input/potential/WallSpring_cavity_radial/equil_dist'][0])
        elif adv.get('cavity_radius'):
            radius = np.full(len(ids), float(adv['cavity_radius']))
        check_spring(pot, 'cavity_radial', 'WallSpring', len(ids), radius, np.full(len(ids), 5.), np.ones(len(ids)))

    if adv.get('ask_before_using_AFM'):
        # coord_const.cpp:59-70, 93: tip k is at start_pos[k] + velocities[k]*t, its start the table's tip position or,
        # for a table without one, its CA in the config
        head, rows = read_table(adv['ask_before_using_AFM'])
        res = [int(r[0]) for r in rows]
        g = pot['MovingConst3D']
        start = (f['input/pos'][:, :, 0][CA[res]] if len(head) == 5
                 else np.array([[float(v) for v in r[2:5]] for r in rows]))
        compare('MovingConst3D', 'id', np.arange(len(rows)), g['id'][:])
        compare('MovingConst3D', 'velocities', np.array([[float(v) for v in r[-3:]] for r in rows]), g['velocities'][:], 1e-12)
        compare('MovingConst3D', 'start_pos', start, g['start_pos'][:], 1e-6)
        compare('Distance3D_moving', 'id (CA, tip row)', np.column_stack((CA[res], np.arange(len(rows)))),
                pot['Distance3D_moving/id'][:].astype(int))
        check_spring(pot, 'pulling', 'Spring', len(rows), np.zeros(len(rows)), [float(r[1]) for r in rows])
    if adv.get('tension'):
        # tension.cpp:26-27, 38-40: row i pulls pos row atom[i] along tension_coeff[i]
        rows = read_table(adv['tension'])[1]
        compare('tension', 'atom (CA of each residue)', CA[[int(r[0]) for r in rows]], pot['tension/atom'][:])
        compare('tension', 'tension_coeff', np.array([[float(v) for v in r[1:4]] for r in rows]),
                pot['tension/tension_coeff'][:], 1e-12)
    if adv.get('contact_energies'):
        # sidechain_radial.cpp:260-263, 273: id rows are rows of placement_fixed_point_only_CB (one per residue)
        rows = read_table(adv['contact_energies'])[1]
        cb_res = pot['placement_fixed_point_only_CB/affine_residue'][:]
        row_of = dict((int(r), i) for i, r in enumerate(cb_res))
        compare('contact', 'id (CB rows of residue1, residue2)', np.array([[row_of[int(r[0])], row_of[int(r[1])]]
                                                                         for r in rows]), pot['contact/id'][:])
        for d, c in (('energy', 2), ('distance', 3), ('width', 4)):
            compare('contact', d, np.array([float(r[c]) for r in rows]), pot['contact/' + d][:], 1e-12)

# ------------------------------------------------------------------------------------------------------------------
# hybrid: layout, dry-MARTINI topology, interface maps
# ------------------------------------------------------------------------------------------------------------------

def read_itp(paths):
    # a minimal GROMACS reader: #define macros, #ifdef/#ifndef/#else/#endif with only FLEXIBLE defined (Upside takes
    # a topology's constraints in their FLEXIBLE harmonic form, up.md "MARTINI-specific nodes"), the sections the
    # checks need
    defines, top = {'FLEXIBLE': []}, dict(atomtypes={}, nonbond={}, molecules={})
    for path in paths:
        live, mol, section = [True], None, None
        for raw in open(path):
            line = raw.split(';')[0].strip()
            if not line: continue
            if line.startswith('#'):
                tok = line.split()
                if tok[0] == '#ifdef': live.append(live[-1] and tok[1] in defines)
                elif tok[0] == '#ifndef': live.append(live[-1] and tok[1] not in defines)
                elif tok[0] == '#else': live[-1] = live[-2] and not live[-1]
                elif tok[0] == '#endif': live.pop()
                elif tok[0] == '#define' and live[-1]: defines[tok[1]] = tok[2:]
                continue
            if not live[-1]: continue
            m = re.match(r'\[\s*(\w+)\s*\]', line)
            if m:
                section = m.group(1).lower(); continue
            tok = sum([defines.get(x, [x]) for x in line.split()], [])
            if section == 'atomtypes':
                top['atomtypes'][tok[0]] = float(tok[1])
            elif section == 'nonbond_params':
                top['nonbond'][(tok[0], tok[1])] = (float(tok[3]), float(tok[4]))   # sigma nm, epsilon kJ/mol
            elif section == 'moleculetype':
                mol = top['molecules'][tok[0].upper()] = dict(nrexcl=int(tok[1]), atoms=[], bonds=[], angles=[],
                                                              dihedrals=[])
            elif section == 'atoms':
                mol['atoms'].append((tok[4].upper(), tok[1], float(tok[6])))
            elif section == 'bonds':
                mol['bonds'].append((int(tok[0])-1, int(tok[1])-1, float(tok[3]), float(tok[4])))
            elif section == 'angles':
                mol['angles'].append((int(tok[0])-1, int(tok[1])-1, int(tok[2])-1, float(tok[4]), float(tok[5])))
            elif section == 'dihedrals':
                # type 1: i j k l funct phi_s k n
                mol['dihedrals'].append(tuple(int(x)-1 for x in tok[:4]) + (float(tok[5]), float(tok[6]), int(tok[7])))
    return top

def env_molecules(inp, start):
    # an environment molecule is a run of atoms sharing residue id and name
    rid, rname = inp['residue_ids'][:], strings(inp['residue_names'])
    n = len(rid)
    cut = [start] + [i for i in range(start+1, n) if rid[i] != rid[i-1] or rname[i] != rname[i-1]] + [n]
    return [(rname[a].upper(), np.arange(a, b)) for a, b in zip(cut[:-1], cut[1:])]

def check_hybrid(f, L, lib, sc):
    inp, pot = f['input'], f['input/potential']
    n, N = L.n_res, L.n_atom
    names, roles = strings(inp['atom_names']), strings(inp['atom_roles'])
    lay = ['N', 'CA', 'C', 'O']*n + ['BB']*n
    # up.md 5 layout against the file's own atom names; martini_hybrid.cpp:394-397 classifies by atom_roles
    compare('/input', 'atom_names[:5n] (N, CA, C, O per residue, then BB)', np.array(lay), np.array(names[:5*n]))
    compare('/input', 'atom_roles[:5n]', np.array(lay), np.array(roles[:5*n]))
    r = np.arange(n)
    bb = inp['hybrid_bb_map']
    # martini_hybrid.cpp:362-375, 75-109: atom_indices (N, CA, C, O), bb_atom_index the BB slot, weights per slot
    compare('hybrid_bb_map', 'atom_indices', np.column_stack((L.N, L.CA, L.C, L.O)), bb['atom_indices'][:])
    compare('hybrid_bb_map', 'atom_mask', np.ones((n, 4), int), bb['atom_mask'][:])
    compare('hybrid_bb_map', 'bb_atom_index', L.BB, bb['bb_atom_index'][:])
    compare('hybrid_bb_map', 'weights (N, CA, C, O masses 14, 12, 12, 16 over 54, up.md 1.4)',
            np.tile(np.array([14., 12., 12., 16.])/54., (n, 1)), bb['weights'][:], 1e-7)
    if 'infer_H_O' in pot:
        # martini_hybrid.cpp:547-557: residue k's O is the acceptor whose slot 1 is atom_indices[k][2]
        acc_c = pot['infer_H_O/acceptors/id'][:, 1]
        compare('hybrid_bb_map', 'residue has an acceptor on its C (all but chain ends)', ~L.last, np.isin(L.C, acc_c))
    memb = inp['hybrid_env_topology/protein_membership'][:]
    compare('hybrid_env_topology', 'protein_membership (one protein, then environment)',
            np.r_[np.zeros(5*n, int), -np.ones(N-5*n, int)], memb)
    compare('hybrid_env_topology', 'env_atom_indices', np.arange(5*n, N), inp['hybrid_env_topology/env_atom_indices'][:])

    top = read_itp([ITP_DIR + x for x in sorted(os.listdir(ITP_DIR)) if x.endswith('.itp')])
    mols = env_molecules(inp, 5*n)
    atype, charge, mass = strings(inp['type']), inp['charges'][:], inp['mass'][:]
    pclass, molid = strings(inp['particle_class']), inp['molecule_ids'][:]
    exp_type, exp_q, exp_m, exp_cls = list(atype[:5*n]), list(charge[:5*n]), list(mass[:5*n]), ['PROTEIN']*(5*n)
    bonds, angles, dihedrals, mol_of, exp_names, got_names = [], [], [], np.zeros(N, int), [], []
    for k, (mname, atoms) in enumerate(mols):
        mt = top['molecules'][mname]
        by_name = dict((a[0], a) for a in mt['atoms'])
        here = [names[a].upper() for a in atoms]
        exp_names.append(mname + ':' + ','.join(sorted(by_name)))
        got_names.append(mname + ':' + ','.join(sorted(here)))
        idx = dict(zip(here, atoms))
        for a, nm in zip(atoms, here):
            t, q = by_name[nm][1:] if nm in by_name else ('?', np.nan)
            exp_type.append(t); exp_q.append(q); exp_m.append(top['atomtypes'].get(t, np.nan)/12.)
            exp_cls.append('LIPID' if len(mt['atoms']) > 1 else 'ION')   # the single-bead molecules are the ions
        mol_of[atoms] = k + 1
        if all(a[0] in idx for a in mt['atoms']):
            at = [idx[a[0]] for a in mt['atoms']]
            for i, j, b0, kb in mt['bonds']:
                bonds.append((at[i], at[j], 10.*b0, kb/E_UP_KJ/100.))       # nm -> A, kJ/mol/nm^2 -> E_up/A^2
            for i, j, k3, th, kt in mt['angles']:
                angles.append((at[i], at[j], at[k3], th, kt/E_UP_KJ))
            for i, j, k3, l, phi_s, kd, n in mt['dihedrals']:
                dihedrals.append((at[i], at[j], at[k3], at[l], phi_s, kd/E_UP_KJ, n))
            if mt['nrexcl'] != 1:
                note('%s has nrexcl %d; the pair check excludes bonded neighbours only' % (mname, mt['nrexcl']))
    compare('/input', 'bead names of each environment molecule (= its itp moleculetype)', np.array(exp_names),
            np.array(got_names))
    compare('/input', 'type (itp bead type; protein part as stored)', np.array(exp_type), np.array(atype))
    compare('/input', 'charges (itp charge)', np.array(exp_q), charge, 1e-6)
    compare('/input', 'mass (itp atomtype mass / 12)', np.array(exp_m), mass, 1e-6)
    compare('/input', 'particle_class', np.array(exp_cls), np.array(pclass))
    cls = np.array(exp_cls)
    compare('/input', 'environment order: lipids, then ions (up.md 5)', True,
            bool(not (cls == 'ION').any() or not (cls == 'LIPID').any() or
                 np.nonzero(cls == 'LIPID')[0].max() < np.nonzero(cls == 'ION')[0].min()))
    # box.cpp:57-73: one barostat molecule per id; protein one molecule, each environment molecule its own
    compare('/input', 'molecule_ids partition (protein, then one id per environment molecule)',
            canonical_labels(mol_of), canonical_labels(molid))
    bonds, angles = np.array(bonds).reshape(-1, 4), np.array(angles).reshape(-1, 5)
    dihedrals = np.array(dihedrals).reshape(-1, 7)

    if 'dist_spring' in pot:
        # martini_potential.cpp:1052-1054, 1100-1110: id (i, j), equil_dist A, spring_const E_up/A^2
        g = pot['dist_spring']
        compare('dist_spring', 'id (sorted pair)', np.sort(bonds[:, :2].astype(int), axis=1), np.sort(g['id'][:], axis=1)
                if g['id'].shape == (len(bonds), 2) else g['id'][:])
        compare('dist_spring', 'equil_dist', bonds[:, 2], g['equil_dist'][:], 1e-5)
        compare('dist_spring', 'spring_const', bonds[:, 3], g['spring_const'][:], 1e-5*max(1., bonds[:, 3].max()))
    if 'angle_spring' in pot:
        # martini_potential.cpp:1158-1160, 1214-1215: the angle is at the MIDDLE atom (GROMACS order i-j-k)
        g = pot['angle_spring']
        e = angles[:, :3].astype(int)
        canon = lambda a: np.column_stack((np.minimum(a[:, 0], a[:, 2]), a[:, 1], np.maximum(a[:, 0], a[:, 2])))
        compare('angle_spring', 'id (vertex in the middle, ends sorted)', canon(e), canon(g['id'][:]) if g['id'].shape == e.shape
                else g['id'][:])
        compare('angle_spring', 'equil_angle_deg', angles[:, 3], g['equil_angle_deg'][:], 1e-4)
        compare('angle_spring', 'spring_const', angles[:, 4], g['spring_const'][:], 1e-6*max(1., angles[:, 4].max()))
    # martini_potential.cpp DihedralSpring: id (i, j, k, l) in the ITP's order, phi_s deg, k E_up, n; written only
    # when some environment molecule's ITP has [ dihedrals ]
    compare('dihedral_spring', 'present exactly when the ITPs give dihedrals', len(dihedrals) > 0, 'dihedral_spring' in pot)
    if 'dihedral_spring' in pot and len(dihedrals):
        g = pot['dihedral_spring']
        compare('dihedral_spring', 'id (ITP order)', dihedrals[:, :4].astype(int), g['id'][:])
        compare('dihedral_spring', 'equil_angle_deg', dihedrals[:, 4], g['equil_angle_deg'][:], 1e-4)
        compare('dihedral_spring', 'spring_const', dihedrals[:, 5], g['spring_const'][:], 1e-6*max(1., dihedrals[:, 5].max()))
        compare('dihedral_spring', 'multiplicity', dihedrals[:, 6].astype(int), g['multiplicity'][:])
    if 'martini_potential' in pot:
        check_martini_pairs(pot['martini_potential'], L, top, bonds, atype, charge)
    if 'martini_sc_table_1body' in pot and sc is not None:
        check_sc_table(pot, L, lib, sc, memb, atype, charge)
    if 'brownian' in inp:
        check_brownian(inp, pot, L, mols)

def canonical_labels(x):
    # relabel ids by order of first appearance, so two partitions compare as arrays
    first = {}
    return np.array([first.setdefault(v, len(first)) for v in x.tolist()])

def check_martini_pairs(g, L, top, bonds, atype, charge):
    n, N = L.n_res, L.n_atom
    # martini_potential.cpp:268-270: atom_indices and charges are read for their size only; the pair charges are the
    # coefficients' (:384)
    compare('martini_potential', 'atom_indices', np.arange(N), g['atom_indices'][:])
    compare('martini_potential', 'charges (= /input/charges)', charge, g['charges'][:], 1e-6)
    # Sites: the BB beads and the environment (up.md 1.3; the engine skips a protein atom other than a derived BB
    # against the environment, :528-537, and every intra-protein pair, :505); bonded neighbours excluded (nrexcl 1)
    sites = np.r_[L.BB, np.arange(5*n, N)].astype(np.int64)
    pairs = g['pairs'][:].astype(np.int64)
    key = pairs[:, 0]*N + pairs[:, 1]
    # :76-79, 107-119, 229-234: first < second, and pairs sorted so each atom's neighbours ascend
    compare('martini_potential', 'pairs first < second', True, bool((pairs[:, 0] < pairs[:, 1]).all()))
    compare('martini_potential', 'pairs sorted (binary-searched neighbour rows)', True, bool((np.diff(key) > 0).all()))
    i, j = np.triu_indices(len(sites), 1)
    expect = sites[i]*N + sites[j]
    del i, j
    b = np.sort(bonds[:, :2].astype(np.int64), axis=1)
    expect = np.setdiff1d(expect, b[:, 0]*N + b[:, 1])
    missing, extra = np.setdiff1d(expect, key), np.setdiff1d(key, expect)
    for k in missing[:MAX_ROWS]: print('MISMATCH martini_potential/pairs missing (%d, %d)' % divmod(k, N))
    for k in extra[:MAX_ROWS]: print('MISMATCH martini_potential/pairs extra (%d, %d)' % divmod(k, N))
    results.append(('martini_potential', 'pairs (BB and environment, bonded excluded)', len(expect), len(missing)+len(extra)))
    del expect
    # coefficients[coefficient_indices[p]] = (eps, sigma, q_first, q_second) of the pair's atoms (:309-329, 381-401)
    cls_key = ['%s|%r' % (t, q) for t, q in zip(atype, charge.tolist())]
    cls_names, cls = np.unique(cls_key, return_inverse=True)
    n_cls = len(cls_names)
    ci = g['coefficient_indices'][:]
    combo = (cls[pairs[:, 0]]*n_cls + cls[pairs[:, 1]])*(ci.max()+1) + ci
    u, first = np.unique(combo, return_index=True)
    coeff = g['coefficients'][:]
    exp_rows, got_rows = [], []
    for c, p in zip(u, first):
        a, bb = pairs[p]
        sig, eps = top['nonbond'].get((atype[a], atype[bb])) or top['nonbond'][(atype[bb], atype[a])]
        exp_rows.append((np.float32(eps/E_UP_KJ), np.float32(10.*sig), charge[a], charge[bb]))
        got_rows.append(coeff[ci[p]])
    compare('martini_potential', 'coefficients[coefficient_indices] (one row per atom-class pair and index)',
            np.array(exp_rows), np.array(got_rows), 1e-6)
    # each class pair must use one coefficient row
    cp = u // (ci.max()+1)
    compare('martini_potential', 'one coefficient row per atom-class pair', len(np.unique(cp)), len(cp))
    # :375-399: each coefficient row's (eps, sigma, q1*q2) is a row of the grid table, exactly in float32
    trip = set(zip(g['unique_eps_eup'][:].astype('f4').tolist(), g['unique_sig_ang'][:].astype('f4').tolist(),
                   g['unique_charge_product'][:].astype('f4').tolist()))
    exact = [(float(c[0]), float(c[1]), float(np.float32(c[2]*c[3]))) in trip for c in coeff]
    compare('martini_potential', 'coefficient row has an exact (eps, sig, qq) grid (no tolerance fallback)',
            np.ones(len(coeff), bool), np.array(exact))

def check_sc_table(pot, L, lib, sc, memb, atype, charge):
    g = pot['martini_sc_table_1body']
    with h5py.File(lib['martini'], 'r') as m:
        sct = m['sc_table']
        restypes = strings(sct['restype_order'])
        ttype, tq = strings(sct['target_type']), sct['target_charge'][:]
        for d in ('grid_ang', 'theta_grid', 'rotamer_full_energy_eup', 'rotamer_count'):
            compare('martini_sc_table_1body', d + ' (= martini.h5 sc_table)', sct[d][:], g[d][:], 0.)
    with tb.open_file(lib['sidechain']) as t:
        names = [x.decode() for x in t.root.restype_order[:]]
        ss = t.root.rotamer_start_stop_bead[:]
        lib_nrot = np.array([(ss[names.index(x)][1]-ss[names.index(x)][0])//ss[names.index(x)][2] for x in restypes])
    compare('martini_sc_table_1body', 'rotamer_count = sidechain library rotamers (martini_build_tables.py _sidechain_conformers)',
            lib_nrot, g['rotamer_count'][:])
    res, rot = sc['res'], sc['rot']
    # martini_potential.cpp:780-785, 898-901: one row per rotamer row of the rotamer node; row_residue_index is the CB
    # row (= residue), row_rotamer_index the table's rotamer (the library's order), row_residue_table_index the
    # table's residue type or -1 for a residue without a MARTINI side chain (:673, :832)
    compare('martini_sc_table_1body', 'row_residue_index', res, g['row_residue_index'][:])
    compare('martini_sc_table_1body', 'row_rotamer_index', rot, g['row_rotamer_index'][:])
    compare('martini_sc_table_1body', 'row_residue_table_index', np.array([restypes.index(L.seq[i]) if L.seq[i] in restypes
                                                                        else -1 for i in res]), g['row_residue_table_index'][:])
    env = np.nonzero(memb < 0)[0]
    # :783-784, 899, 904: env_atom_index the atom, env_target_index the table target of its (bead type, charge)
    compare('martini_sc_table_1body', 'env_atom_index (every environment atom)', env, g['env_atom_index'][:])
    tgt = dict(((t, float(q)), k) for k, (t, q) in enumerate(zip(ttype, tq)))
    compare('martini_sc_table_1body', 'env_target_index (bead type and charge)',
            np.array([tgt.get((atype[a], float(charge[a])), -1) for a in env]), g['env_target_index'][:])

def check_brownian(inp, pot, L, mols):
    br = inp['brownian']
    n = L.n_res
    carriers = np.sort(np.r_[L.N, L.CA, L.C])
    env = np.arange(5*n, L.n_atom)
    # martini_brownian.cpp:52-88: atom_index lists the atoms the Langevin step drives, friction[i] belongs to
    # atom_index[i]; up.md 1.4: every dynamic particle, the derived O and BB slots absent
    atom_index = br['atom_index'][:]
    compare('brownian', 'atom_index (N, CA, C and the environment)', np.r_[carriers, env], atom_index)
    phase = br.attrs['dynamics_phase'].decode()
    if phase != 'production' or len(atom_index) != len(carriers) + len(env):
        note('brownian friction not checked (phase %s)' % phase); return
    fr = br['friction'][:]
    lipid_f = dict(zip([x.decode() for x in br.attrs['lipid_names']], np.atleast_1d(br.attrs['lipid_bead_friction_up'])))
    bead_f, is_lipid = np.zeros(L.n_atom), np.zeros(L.n_atom, bool)
    for mname, atoms in mols:
        bead_f[atoms] = lipid_f[mname] if mname in lipid_f else br.attrs['ion_bead_friction_up']
        is_lipid[atoms] = mname in lipid_f
    compare('brownian', 'friction of environment atoms (their lipid, or the ion value)', bead_f[env],
            fr[len(carriers):], 1e-6)
    # A carrier's friction is the sum of the frictions of the lipid beads within interface_friction_cutoff_A, counted
    # at a stage handoff (martini_prepare_system.py:1220-1274); protein_contact_atom_index lists the carriers with a
    # contact, protein_lipid_contact_count their counts.  A continuation keeps its handoff's friction.
    contact, count = br['protein_contact_atom_index'][:], br['protein_lipid_contact_count'][:]
    compare('brownian', 'protein_contact_atom_index are carriers', np.ones(len(contact), bool), np.isin(contact, carriers))
    stored_count = np.zeros(L.n_atom, int); stored_count[contact] = count
    if len(lipid_f) == 1:
        compare('brownian', 'friction of carriers = contact count x lipid bead friction',
                stored_count[carriers]*list(lipid_f.values())[0], fr[:len(carriers)], 1e-4)
    x = inp['pos'][:, :, 0].astype('f8')
    box = np.array([pot['martini_potential'].attrs[k] for k in ('x_len', 'y_len', 'z_len')], dtype='f8')
    d = x[carriers][:, None, :] - x[is_lipid][None, :, :]
    d -= box*np.round(d/box)
    within = (d**2).sum(axis=-1) < br.attrs['interface_friction_cutoff_A']**2
    agree = within.sum(axis=1) == stored_count[carriers]
    print('INFO brownian carrier contacts recounted at /input/pos agree for %d of %d carriers' % (agree.sum(), len(agree)))
    if agree.all():
        compare('brownian', 'friction of carriers (recounted at /input/pos: written at these positions)',
                (within*bead_f[is_lipid][None, :]).sum(axis=1), fr[:len(carriers)], 1e-4)
    else:
        note('brownian carrier counts are not those of /input/pos (%d of %d differ): the friction was counted at an '
             'earlier handoff' % ((~agree).sum(), len(agree)))

# ------------------------------------------------------------------------------------------------------------------

NO_INDEX = {'martini_hybrid_position': 'no datasets; its maps are /input/hybrid_bb_map'}

def main():
    for path in sys.argv[1:]:
        del results[:]; del notes[:]
        print('==== %s' % path)
        args = read_args(path)
        with h5py.File(path, 'r') as f:
            L = Layout(f)
            lib = libraries(args, L.hybrid)
            pot = f['input/potential']
            print('INFO %d residues, %d atoms, stride %d, chain starts %s' % (L.n_res, L.n_atom, L.stride,
                                                                            list(L.chain_starts)))
            check_backbone(pot, L, args)
            check_rama(pot, L, lib, args)
            check_pivots(f, pot, L)
            check_hbond(pot, L, args, lib)
            sc = check_sidechain(pot, L, args, lib)
            check_cb(pot, L, args)
            check_environment(pot, L, args, lib, sc)
            check_membrane(pot, L, args, lib)
            check_advanced(f, pot, L, read_args(path, 'adv_args'))
            if L.hybrid: check_hybrid(f, L, lib, sc)
            seen = set(x[0] for x in results)
            unchecked = [g for g in pot if g not in seen]
        n_bad = sum(1 for x in results if x[3] != 0)
        print('---- summary %s' % os.path.basename(path))
        for node, dset, n_rows, bad in results:
            print('%-4s %-40s %-90s %8s rows  %s' % ('ok' if bad == 0 else 'BAD', node, dset, n_rows,
                                                    '' if bad == 0 else '%s differ' % bad))
        for g in unchecked:
            print('not checked: %s (%s)' % (g, NO_INDEX.get(g, 'no check written')))
        print('==== %s: %d of %d datasets mismatch' % (os.path.basename(path), n_bad, len(results)))

if __name__ == '__main__':
    main()
