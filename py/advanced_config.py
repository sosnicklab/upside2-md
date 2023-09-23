#!/usr/bin/env python

from itertools import count, groupby
import numpy as np
import tables as tb
import sys,os
from upside_nodes import *

three_letter_aa = dict(
        A='ALA', C='CYS', D='ASP', E='GLU',
        F='PHE', G='GLY', H='HIS', I='ILE',
        K='LYS', L='LEU', M='MET', N='ASN',
        P='PRO', Q='GLN', R='ARG', S='SER',
        T='THR', V='VAL', W='TRP', Y='TYR')

aa_num = dict([(k,i) for i,k in enumerate(sorted(three_letter_aa.values()))])
one_letter_aa = dict([(v,k) for k,v in three_letter_aa.items()])
one_letter_aa['CPR'] = 'P' # add entry for cis proline
deg=np.deg2rad(1)
default_filter = tb.Filters(complib='zlib', complevel=5, fletcher32=True)
n_bit_rotamer = 4

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def bstring(string):
    return bytes(string, encoding="ascii")

def highlight_residues(name, fasta, residues_to_highlight):
    fasta_one_letter = [one_letter_aa[x.decode('ASCII')] for x in fasta]
    residues_to_highlight = set(residues_to_highlight)
    print ('%s:  %s' % (name, ''.join((f.upper() if i in residues_to_highlight else f.lower()) for i,f in enumerate(fasta_one_letter))))

def vmag(x):
    assert x.shape[-1] == 3
    return np.sqrt(x[...,0]**2+x[...,1]**2+x[...,2]**2)


def chain_endpts(n_res, chain_first_residue, i):
    n_chains = chain_first_residue.size+1
    if i == 0:
        first_res = 0
        next_first_res = chain_first_residue[i]
    elif i == n_chains-1:
        first_res = chain_first_residue[i-1]
        next_first_res = n_res
    else:
        first_res = chain_first_residue[i-1]
        next_first_res = chain_first_residue[i]

    return first_res, next_first_res

def create_array(grp, nm, obj=None):
    return t.create_earray(grp, nm, obj=obj, filters=default_filter)

def parse_residue_list(input_string):
    num_list = []
    for ss in input_string.split(','):
        sss = ss.split('-')
        if len(sss) > 1:
            num_list += range(int(sss[0]), int(sss[1])+1)
        else:
            num_list.append( int(sss[0]))
    return num_list

def add_to_const3D1(xyz):
    pot_group = t.root.input.potential

    if 'Const3D' in pot_group:
        xyz  = np.array(xyz)
        xyz0 = pot_group.Const3D.value[:]
        id0  = pot_group.Const3D.id[:]
        index = -1
        for i,x in enumerate(xyz0):
            if np.sum((x-xyz)**2)**0.5 < 1e-6:
                index = i*1
                break
        if index < 0:
            index = xyz0.shape[0]
            new_xyz = np.concatenate((xyz0,  [xyz]))
            new_id  = np.concatenate((id0, [index]))
            pot_group.Const3D.value._f_remove(recursive=True, force=True)
            pot_group.Const3D.id._f_remove(recursive=True, force=True)
            create_array(pot_group.Const3D, 'value', new_xyz)
            create_array(pot_group.Const3D, 'id',    new_id)
    else:
        xyz     = np.array([xyz])
        n_const = xyz.shape[0]
        const   = Const3DCoord1(t, None, 'pos', 0, 1, 2, xyz)
        index   = 0

    return index

def add_to_const3D2(new_id, input_name):
    pot_group = t.root.input.potential

    node_name = 'Const3D_{}'.format(input_name)

    if node_name in pot_group:

        node = pot_group[node_name]

        id0  = node.id[:]
        index = -1
        for i,x in enumerate(id0):
            if x == new_id:
                index = i*1
                break
        if index < 0:
            index = id0.size
            new_ids  = np.concatenate((id0, [new_id]))
            node.id._f_remove(recursive=True, force=True)
            create_array(node, 'id', new_ids)
    else:
        const   = Const3DCoord2(t, input_name, 0, 1, 2, [new_id])
        index   = 0

    return index

def append_const3D_to_pos():
    ConcatCoord(t, 'pos_const', ['pos', 'Const3D'])

def write_affine_alignment(n_res):
    grp = t.create_group(potential, 'affine_alignment')
    grp._v_attrs.arguments = np.array([b'pos'])

    ref_geom = np.zeros((n_res,3,3))
    ref_geom[:,0] = (-1.19280531, -0.83127186, 0.)  # N
    ref_geom[:,1] = ( 0.,          0.,         0.)  # CA
    ref_geom[:,2] = ( 1.25222632, -0.87268266, 0.)  # C
    ref_geom -= ref_geom.mean(axis=1)[:,None]

    N  = np.arange(n_res)*3 + 0
    CA = np.arange(n_res)*3 + 1
    C  = np.arange(n_res)*3 + 2

    atoms = np.column_stack((N,CA,C))
    create_array(grp, 'atoms', obj=atoms)
    create_array(grp, 'ref_geom', obj=ref_geom)

def write_CB(fasta):
    # Place CB
    pgrp = t.create_group(potential, 'placement_fixed_point_only_CB')
    pgrp._v_attrs.arguments = np.array([b'affine_alignment'])
    ref_pos = np.zeros((4,3))
    ref_pos[0] = (-1.19280531, -0.83127186,  0.)        # N
    ref_pos[1] = ( 0.,          0.,          0.)        # CA
    ref_pos[2] = ( 1.25222632, -0.87268266,  0.)        # C
    ref_pos[3] = ( 0.,          0.94375626,  1.2068012) # CB
    ref_pos -= ref_pos[:3].mean(axis=0,keepdims=1)

    placement_data = np.zeros((1,3))
    placement_data[0,0:3] = ref_pos[3]

    create_array(pgrp, 'affine_residue',  np.arange(len(fasta)))
    create_array(pgrp, 'layer_index',     np.zeros(len(fasta),dtype='i'))
    create_array(pgrp, 'placement_data',  placement_data)

#=========================================================================
#                             wall potential 
#=========================================================================

def write_cavity_radial(radius, spring_const=5.):
    # add a (0,0,0) point
    const_index  = add_to_const3D1([0, 0, 0])
    # distance between pos and (0,0,0)
    ids = []
    for i in range(n_atom):
        ids.append([i, const_index])
    ids = np.array(ids)
    cavity_radial_dist = DistanceCoord(t, 'cavity', 'pos', 'Const3D', ids )
    # make wall spring potential
    ids = np.arange(n_atom)
    R0  = np.ones(n_atom)*radius
    K   = np.ones(n_atom)*spring_const
    W   = np.ones(n_atom) # 0 means left wall, 1 means right wall
    cavity_radial_node = WallSpringPotential(t, 'cavity_radial', cavity_radial_dist, 0, ids, R0, K, W)

def write_const_wall_spring(parser, fasta, spring_table):

    header_basic = 'residue radius spring_const wall_type'

    fields        = [ ln.split() for ln in open(spring_table) ]
    actual_header = [  x.lower() for  x in fields[0] ]

    z_only = False
    xy_only = False
    all_xyz = False
    if actual_header   == (header_basic + ' x0').split():
        dim = 0; z_only = True
    elif actual_header == (header_basic + ' y0').split():
        dim = 1; z_only = True
    elif actual_header == (header_basic + ' z0').split():
        dim = 2; z_only = True
    elif actual_header == (header_basic + ' x0 y0').split():
        dim1 = 0; dim2 = 1; xy_only = True
    elif actual_header == (header_basic + ' x0 z0').split():
        dim1 = 0; dim2 = 2; xy_only = True
    elif actual_header == (header_basic + ' y0 z0').split():
        dim1 = 1; dim2 = 2; xy_only = True
    elif actual_header == (header_basic + ' x0 y0 z0').split():
        all_xyz = True
    else:
        parser.error('Format error of the first line')

    if not all(len(f)==len(fields[0]) for f in fields):
        parser.error('Invalid format for %s file' % spring_table)

    fields = fields[1:]

    atom         = []
    radius       = []
    spring_const = []
    wall_type    = []
    xyz          = []
    for i,f in enumerate(fields):
        group = parse_residue_list(f[0])
        atom += group

        for j in range(len(group)):
            msg = '%s energy specified for residue %i (zero is first residue) but there are only %i residues in the FASTA'
            if not (0 <= group[j] < len(fasta)): raise ValueError(msg % (spring_table, group[j], len(fasta)))
            radius.append(float(f[1]))
            spring_const.append(float(f[2]))
            wall_type.append(float(f[3]))
            xxx = [0., 0., 0.]
            if z_only:
                xxx[dim]  = float(f[4])
            elif xy_only:
                xxx[dim1] = float(f[4])
                xxx[dim2] = float(f[5])
            else:
                for k in range(3):
                    xxx[k] = float(f[4+k])
            xyz.append(xxx)
    atom         = np.array(atom)
    atom         = atom*3+1  # Ca atoms
    radius       = np.array(radius)
    spring_const = np.array(spring_const)
    wall_type    = np.array(wall_type)
    xyz          = np.array(xyz)

    n_spring     = atom.shape[0]

    ids = []
    for i in range(n_spring):
        ndx = add_to_const3D1(xyz[i])
        ids.append([atom[i], ndx])
    ids = np.array(ids)

    pot_group = t.root.input.potential
    if z_only:
        flat_bottom_dist = Distance1DCoord(t, 'wall', 'pos', 'Const3D', dim, ids)
    elif xy_only:
        flat_bottom_dist = Distance2DCoord(t, 'wall', 'pos', 'Const3D', dim1, dim2, ids)
    else:
        flat_bottom_dist = DistanceCoord  (t, 'wall', 'pos', 'Const3D', ids)

    # make wall spring potential
    ids = np.arange(n_spring)
    used_dim = 0
    flat_bottom_node = WallSpringPotential(t, 'wall', flat_bottom_dist, used_dim, ids, radius, spring_const, wall_type)

def write_pair_wall_spring(parser, fasta, spring_table):

    header_basic = 'residue1 residue2 radius spring_const wall_type'

    fields        = [ ln.split() for ln in open(spring_table) ]
    actual_header = [  x.lower() for  x in fields[0] ]

    z_only = False
    xy_only = False
    all_xyz = False
    if actual_header   == (header_basic + ' x').split():
        dim = 0; z_only = True
    elif actual_header == (header_basic + ' y').split():
        dim = 1; z_only = True
    elif actual_header == (header_basic + ' z').split():
        dim = 2; z_only = True
    elif actual_header == (header_basic + ' xy').split():
        dim1 = 0; dim2 = 1; xy_only = True
    elif actual_header == (header_basic + ' xz').split():
        dim1 = 0; dim2 = 2; xy_only = True
    elif actual_header == (header_basic + ' yz').split():
        dim1 = 1; dim2 = 2; xy_only = True
    elif actual_header == (header_basic + ' xyz').split():
        all_xyz = True
    else:
        parser.error('Format error of the first line')

    if not all(len(f)==(len(fields[0])-1) for f in fields[1:]):
        parser.error('Invalid format for %s file' % spring_table)

    fields = fields[1:]

    atom         = []
    radius       = []
    spring_const = []
    wall_type    = []
    xyz          = []
    for i,f in enumerate(fields):
        atom.append([int(f[0]), int(f[1])])
        msg = '%s energy specified for residue %i (zero is first residue) but there are only %i residues in the FASTA'
        if not (0 <= atom[-1][0] < len(fasta)): raise ValueError(msg % (spring_table, atom[-1][0], len(fasta)))
        if not (0 <= atom[-1][1] < len(fasta)): raise ValueError(msg % (spring_table, atom[-1][1], len(fasta)))
        radius.append(float(f[2]))
        spring_const.append(float(f[3]))
        wall_type.append(float(f[4]))
    atom         = np.array(atom)
    atom         = atom*3+1  # Ca atoms
    radius       = np.array(radius)
    spring_const = np.array(spring_const)
    wall_type    = np.array(wall_type)
    n_spring     = atom.shape[0]

    pot_group = t.root.input.potential
    if z_only:
        flat_bottom_dist = Distance1DCoord(t, 'pair_wall', 'pos', 'pos', dim, atom)
    elif xy_only:
        flat_bottom_dist = Distance2DCoord(t, 'pair_wall', 'pos', 'pos', dim1, dim2, atom)
    else:
        flat_bottom_dist = DistanceCoord  (t, 'pair_wall', 'pos', 'pos', atom)

    # make wall spring potential
    ids = np.arange(n_spring)
    used_dim = 0
    flat_bottom_node = WallSpringPotential(t, 'pair_wall', flat_bottom_dist, used_dim, ids, radius, spring_const, wall_type)

#=========================================================================
#                            spring potential 
#=========================================================================

def write_const_spring(parser, fasta, spring_table):

    header_basic = 'residue radius spring_const'

    fields        = [ ln.split() for ln in open(spring_table) ]
    actual_header = [  x.lower() for  x in fields[0] ]

    z_only = False
    xy_only = False
    all_xyz = False
    if actual_header   == (header_basic + ' x0').split():
        dim = 0; z_only = True
    elif actual_header == (header_basic + ' y0').split():
        dim = 1; z_only = True
    elif actual_header == (header_basic + ' z0').split():
        dim = 2; z_only = True
    elif actual_header == (header_basic + ' x0 y0').split():
        dim1 = 0; dim2 = 1; xy_only = True
    elif actual_header == (header_basic + ' x0 z0').split():
        dim1 = 0; dim2 = 2; xy_only = True
    elif actual_header == (header_basic + ' y0 z0').split():
        dim1 = 1; dim2 = 2; xy_only = True
    elif actual_header == (header_basic + ' x0 y0 z0').split():
        all_xyz = True
    else:
        parser.error('Format error of the first line')

    if not all(len(f)==len(fields[0]) for f in fields):
        parser.error('Invalid format for %s file' % spring_table)

    fields = fields[1:]

    atom         = []
    radius       = []
    spring_const = []
    xyz          = []
    for i,f in enumerate(fields):
        group = parse_residue_list(f[0])
        atom += group

        for j in range(len(group)):
            msg = '%s energy specified for residue %i (zero is first residue) but there are only %i residues in the FASTA'
            if not (0 <= group[j] < len(fasta)): raise ValueError(msg % (spring_table, group[j], len(fasta)))
            radius.append(float(f[1]))
            spring_const.append(float(f[2]))
            xxx = [0., 0., 0.]
            if z_only:
                xxx[dim]  = float(f[3])
            elif xy_only:
                xxx[dim1] = float(f[3])
                xxx[dim2] = float(f[4])
            else:
                for k in range(3):
                    xxx[k] = float(f[3+k])
            xyz.append(xxx)
    atom         = np.array(atom)
    atom         = atom*3+1  # Ca atoms
    radius       = np.array(radius)
    spring_const = np.array(spring_const)
    xyz          = np.array(xyz)

    n_spring     = atom.shape[0]

    ids = []
    for i in range(n_spring):
        ndx = add_to_const3D1(xyz[i])
        ids.append([atom[i], ndx])
    ids = np.array(ids)

    pot_group = t.root.input.potential
    if z_only:
        spring_dist = Distance1DCoord(t, 'spring', 'pos', 'Const3D', dim, ids)
    elif xy_only:
        spring_dist = Distance2DCoord(t, 'spring', 'pos', 'Const3D', dim1, dim2, ids)
    else:
        spring_dist = DistanceCoord  (t, 'spring', 'pos', 'Const3D', ids)

    # make spring potential
    ids = np.arange(n_spring)
    used_dim = 0
    is_pbc = 0
    box_len = 0.0
    spring_node = SpringPotential(t, 'spring', spring_dist, used_dim, is_pbc, box_len, ids, radius, spring_const)

def write_pair_spring(parser, fasta, spring_table):

    header_basic = 'residue1 residue2 radius spring_const'

    fields        = [ ln.split() for ln in open(spring_table) ]
    actual_header = [  x.lower() for  x in fields[0] ]

    z_only = False
    xy_only = False
    all_xyz = False
    if actual_header   == (header_basic + ' x').split():
        dim = 0; z_only = True
    elif actual_header == (header_basic + ' y').split():
        dim = 1; z_only = True
    elif actual_header == (header_basic + ' z').split():
        dim = 2; z_only = True
    elif actual_header == (header_basic + ' xy').split():
        dim1 = 0; dim2 = 1; xy_only = True
    elif actual_header == (header_basic + ' xz').split():
        dim1 = 0; dim2 = 2; xy_only = True
    elif actual_header == (header_basic + ' yz').split():
        dim1 = 1; dim2 = 2; xy_only = True
    elif actual_header == (header_basic + ' xyz').split():
        all_xyz = True
    else:
        parser.error('Format error of the first line')

    if not all(len(f)==(len(fields[0])-1) for f in fields[1:]):
        parser.error('Invalid format for %s file' % spring_table)

    fields = fields[1:]

    atom         = []
    radius       = []
    spring_const = []
    xyz          = []
    for i,f in enumerate(fields):
        atom.append([int(f[0]), int(f[1])])
        msg = '%s energy specified for residue %i (zero is first residue) but there are only %i residues in the FASTA'
        if not (0 <= atom[-1][0] < len(fasta)): raise ValueError(msg % (spring_table, atom[-1][0], len(fasta)))
        if not (0 <= atom[-1][1] < len(fasta)): raise ValueError(msg % (spring_table, atom[-1][1], len(fasta)))
        radius.append(float(f[2]))
        spring_const.append(float(f[3]))
    atom         = np.array(atom)
    atom         = atom*3+1  # Ca atoms
    radius       = np.array(radius)
    spring_const = np.array(spring_const)

    n_spring     = atom.shape[0]

    if z_only:
        spring_dist = Distance1DCoord(t, 'pair_spring', 'pos', 'pos', dim, atom)
    elif xy_only:
        spring_dist = Distance2DCoord(t, 'pair_spring', 'pos', 'pos', dim1, dim2, atom)
    else:
        spring_dist = DistanceCoord  (t, 'pair_spring', 'pos', 'pos', atom)

    # make wall spring potential
    ids = np.arange(n_spring)
    used_dim = 0
    is_pbc = 0
    box_len = 0.0
    spring_node = SpringPotential(t, 'pair_spring', spring_dist, used_dim, is_pbc, box_len, ids, radius, spring_const)

def write_nail_spring(parser, fasta, spring_table):

    header_basic = 'residue spring_const'

    fields        = [ ln.split() for ln in open(spring_table) ]
    actual_header = [  x.lower() for  x in fields[0] ]

    z_only  = False
    xy_only = False
    all_xyz = False
    if actual_header   == (header_basic + ' x').split():
        dim = 0; z_only = True
    elif actual_header == (header_basic + ' y').split():
        dim = 1; z_only = True
    elif actual_header == (header_basic + ' z').split():
        dim = 2; z_only = True
    elif actual_header == (header_basic + ' xy').split():
        dim1 = 0; dim2 = 1; xy_only = True
    elif actual_header == (header_basic + ' xz').split():
        dim1 = 0; dim2 = 2; xy_only = True
    elif actual_header == (header_basic + ' yz').split():
        dim1 = 1; dim2 = 2; xy_only = True
    elif actual_header == (header_basic + ' xyz').split():
        all_xyz = True
    else:
        parser.error('Format error of the first line')

    if not all(len(f)==(len(fields[0])-1) for f in fields[1:]):
        parser.error('Invalid format for %s file' % spring_table)

    fields = fields[1:]

    atom         = []
    radius       = []
    spring_const = []
    xyz          = []
    for i,f in enumerate(fields):
        group = parse_residue_list(f[0])
        atom += group

        for j in range(len(group)):
            msg = '%s energy specified for residue %i (zero is first residue) but there are only %i residues in the FASTA'
            if not (0 <= group[j] < len(fasta)): raise ValueError(msg % (spring_table, group[j], len(fasta)))
            radius.append(0.0)
            spring_const.append(float(f[1]))
    atom         = np.array(atom)
    atom         = atom*3+1  # Ca atoms
    radius       = np.array(radius)
    spring_const = np.array(spring_const)

    n_spring     = atom.shape[0]

    ids = []
    for i in range(n_spring):
        index = add_to_const3D2(atom[i], 'pos')
        ids.append([atom[i], index])
    ids = np.array(ids)

    if z_only:
        nail_dist = Distance1DCoord(t, 'nail', 'pos', 'Const3D_pos', dim,        ids)
    elif xy_only:
        nail_dist = Distance2DCoord(t, 'nail', 'pos', 'Const3D_pos', dim1, dim2, ids)
    else:
        nail_dist = DistanceCoord  (t, 'nail', 'pos', 'Const3D_pos',             ids)

    # make spring potential
    ids       = np.arange(n_spring)
    used_dim  = 0
    is_pbc = 0
    box_len = 0.0
    nail_node = SpringPotential(t, 'nail', nail_dist, used_dim, is_pbc, box_len, ids, radius, spring_const)

def make_restraint_group(group_num, residues, initial_pos, strength):
    np.random.seed(314159)  # make groups deterministic

    r_atoms = np.array([(3*i+0,3*i+1,3*i+2) for i in sorted(residues)]).reshape((-1,))
    random_pairing = lambda: np.column_stack((r_atoms, np.random.permutation(r_atoms)))

    pairs      = np.concatenate([random_pairing() for i in range(2)], axis=0)
    pairs      = [((x,y) if x<y else (y,x)) for x,y in pairs if x/3!=y/3]   # avoid same-residue restraints
    pairs      = np.array(sorted(set(pairs)))
    pair_dists = vmag(initial_pos[pairs[:,0]]-initial_pos[pairs[:,1]])

    restraint_dist = DistanceCoord  (t, 'restraint{}'.format(group_num), 'pos', 'pos', pairs)
    used_dim = 0
    is_pbc   = 0
    box_len  = 0.0
    ids            = np.arange(len(pairs))
    spring_const   = strength*np.ones(len(pairs))
    restraint_node = SpringPotential(t, 'restraint{}'.format(group_num), restraint_dist, used_dim, is_pbc, box_len, ids, pair_dists, spring_const)


def write_group_node(parser, fasta, spring_table):

    fields        = [ ln.split() for ln in open(spring_table) ]
    actual_header = [  x.lower() for  x in fields[0] ]

    if not all(len(f)==(len(fields[0])) for f in fields):
        parser.error('Invalid format for %s file' % spring_table)

    fields = fields[1:]
    nodes = [0]
    index = []
    for i,f in enumerate(fields):
        atom = parse_residue_list(f[0])
        index += atom
        nodes.append(len(index))

    index = np.array(index)

    msg = '%s specified for residue %i (zero is first residue) but there are only %i residues in the FASTA'
    for j in index:
        if not (0 <= j < len(fasta)): raise ValueError(msg % (spring_table, j, len(fasta)))

    index = index*3+1  # Ca atoms
    nodes = np.array(nodes)

    node = COMCoord(t, 'xyz', 'pos', 3, index, nodes, np.array([0,1,2]))

def write_group_dist_spring(parser, fasta, spring_table):

    header_basic = 'group1 group2 radius spring_const'

    fields        = [ ln.split() for ln in open(spring_table) ]
    actual_header = [  x.lower() for  x in fields[0] ]

    z_only = False
    xy_only = False
    all_xyz = False
    if actual_header   == (header_basic + ' x').split():
        dim = 0; z_only = True
    elif actual_header == (header_basic + ' y').split():
        dim = 1; z_only = True
    elif actual_header == (header_basic + ' z').split():
        dim = 2; z_only = True
    elif actual_header == (header_basic + ' xy').split():
        dim1 = 0; dim2 = 1; xy_only = True
    elif actual_header == (header_basic + ' xz').split():
        dim1 = 0; dim2 = 2; xy_only = True
    elif actual_header == (header_basic + ' yz').split():
        dim1 = 1; dim2 = 2; xy_only = True
    elif actual_header == (header_basic + ' xyz').split():
        all_xyz = True
    else:
        parser.error('Format error of the first line')

    if not all(len(f)==(len(fields[0])-1) for f in fields[1:]):
        parser.error('Invalid format for %s file' % spring_table)

    #assert all(len(f)==(len(fields[0])-1) for f in fields)

    fields = fields[1:]

    atom         = []
    radius       = []
    spring_const = []
    xyz          = []
    for i,f in enumerate(fields):
        atom.append([int(f[0]), int(f[1])])
        radius.append(float(f[2]))
        spring_const.append(float(f[3]))
    atom         = np.array(atom)
    radius       = np.array(radius)
    spring_const = np.array(spring_const)
    n_spring = atom.shape[0]

    if z_only:
        spring_dist = Distance1DCoord(t, 'group_dist', 'GroupCenter_xyz', 'GroupCenter_xyz', dim, atom)
    elif xy_only:
        spring_dist = Distance2DCoord(t, 'group_dist', 'GroupCenter_xyz', 'GroupCenter_xyz', dim1, dim2, atom)
    else:
        spring_dist = DistanceCoord  (t, 'group_dist', 'GroupCenter_xyz', 'GroupCenter_xyz', atom)

    # make wall spring potential
    ids = np.arange(n_spring)
    used_dim = 0
    is_pbc = 0
    box_len = 0.0
    spring_node = SpringPotential(t, 'group_dist', spring_dist, used_dim, is_pbc, box_len, ids, radius, spring_const)

def write_group_const_spring(parser, fasta, spring_table):

    header_basic = 'group radius spring_const'

    fields        = [ ln.split() for ln in open(spring_table) ]
    actual_header = [  x.lower() for  x in fields[0] ]

    z_only = False
    xy_only = False
    all_xyz = False
    if actual_header   == (header_basic + ' x0').split():
        dim = 0; z_only = True
    elif actual_header == (header_basic + ' y0').split():
        dim = 1; z_only = True
    elif actual_header == (header_basic + ' z0').split():
        dim = 2; z_only = True
    elif actual_header == (header_basic + ' x0 y0').split():
        dim1 = 0; dim2 = 1; xy_only = True
    elif actual_header == (header_basic + ' x0 z0').split():
        dim1 = 0; dim2 = 2; xy_only = True
    elif actual_header == (header_basic + ' y0 z0').split():
        dim1 = 1; dim2 = 2; xy_only = True
    elif actual_header == (header_basic + ' x0 y0 z0').split():
        all_xyz = True
    else:
        parser.error('Format error of the first line')

    if not all(len(f)==len(fields[0]) for f in fields):
        parser.error('Invalid format for %s file' % spring_table)

    fields = fields[1:]

    atom         = []
    radius       = []
    spring_const = []
    xyz          = []
    for i,f in enumerate(fields):

        atom.append(float(f[0]))
        radius.append(float(f[1]))
        spring_const.append(float(f[2]))
        xxx = [0., 0., 0.]
        if z_only:
            xxx[dim]  = float(f[3])
        elif xy_only:
            xxx[dim1] = float(f[3])
            xxx[dim2] = float(f[4])
        else:
            for k in range(3):
                xxx[k] = float(f[3+k])
        xyz.append(xxx)
    atom         = np.array(atom)
    radius       = np.array(radius)
    spring_const = np.array(spring_const)
    xyz          = np.array(xyz)

    n_spring     = atom.shape[0]

    ids = []
    for i in range(n_spring):
        ndx = add_to_const3D1(xyz[i])
        ids.append([atom[i], ndx])
    ids = np.array(ids)

    #const3d_group_const = Const3DCoord(t, 'group_const', 'GroupCenter_xyz', 0, 1, 2, np.arange(n_spring), xyz)

    pot_group = t.root.input.potential
    if z_only:
        spring_dist = Distance1DCoord(t, 'group_const', 'GroupCenter_xyz', 'Const3D', dim, ids)
    elif xy_only:
        spring_dist = Distance2DCoord(t, 'group_const', 'GroupCenter_xyz', 'Const3D', dim1, dim2, ids)
    else:
        spring_dist = DistanceCoord  (t, 'group_const', 'GroupCenter_xyz', 'Const3D', ids)

    # make spring potential
    ids = np.arange(n_spring)
    used_dim    = 0
    is_pbc      = 0
    box_len     = 0.0
    spring_node = SpringPotential(t, 'spring', spring_dist, used_dim, is_pbc, box_len, ids, radius, spring_const )

def write_group_nail_spring(parser, fasta, spring_table):

    header_basic = 'group spring_const'

    fields        = [ ln.split() for ln in open(spring_table) ]
    actual_header = [  x.lower() for  x in fields[0] ]

    z_only  = False
    xy_only = False
    all_xyz = False
    if actual_header   == (header_basic + ' x').split():
        dim = 0; z_only = True
    elif actual_header == (header_basic + ' y').split():
        dim = 1; z_only = True
    elif actual_header == (header_basic + ' z').split():
        dim = 2; z_only = True
    elif actual_header == (header_basic + ' xy').split():
        dim1 = 0; dim2 = 1; xy_only = True
    elif actual_header == (header_basic + ' xz').split():
        dim1 = 0; dim2 = 2; xy_only = True
    elif actual_header == (header_basic + ' yz').split():
        dim1 = 1; dim2 = 2; xy_only = True
    elif actual_header == (header_basic + ' xyz').split():
        all_xyz = True
    else:
        parser.error('Format error of the first line')

    if not all(len(f)==(len(fields[0])-1) for f in fields[1:]):
        parser.error('Invalid format for %s file' % spring_table)

    fields = fields[1:]

    atom         = []
    radius       = []
    spring_const = []
    xyz          = []
    for i,f in enumerate(fields):
        atom.append(float(f[0]))
        radius.append(0.0)
        spring_const.append(float(f[1]))
    atom         = np.array(atom)
    radius       = np.array(radius)
    spring_const = np.array(spring_const)

    n_spring     = atom.shape[0]

    const3d_group_nail = Const3DCoord(t, 'GroupCenter_xyz', 0, 1, 2, np.arange(n_spring))

    ids = []
    for i in range(n_spring):
        ids.append([atom[i], i])
    ids = np.array(ids)

    if z_only:
        group_nail_dist = Distance1DCoord(t, 'group_nail', 'GroupCenter_xyz', const3d_group_nail, dim,        ids)
    elif xy_only:
        group_nail_dist = Distance2DCoord(t, 'group_nail', 'GroupCenter_xyz', const3d_group_nail, dim1, dim2, ids)
    else:
        group_nail_dist = DistanceCoord  (t, 'group_nail', 'GroupCenter_xyz', const3d_group_nail,             ids)

    # make spring potential
    ids       = np.arange(n_spring)
    used_dim  = 0
    is_pbc    = 0
    box_len   = 0.0
    nail_node = SpringPotential(t, 'group_nail', group_nail_dist, used_dim, is_pbc, box_len, ids, radius, spring_const)


#=========================================================================
#                          pulling simulation 
#=========================================================================


def write_pulling(parser, fasta, AFM_table, time_initial, time_step):
    fields = [ln.split() for ln in open(AFM_table)]
    header1 = 'residue spring_const tip_pos_x tip_pos_y tip_pos_z pulling_vel_x pulling_vel_y pulling_vel_z'
    header2 = 'residue spring_const pulling_vel_x pulling_vel_y pulling_vel_z'
    actual_header = [x.lower() for x in fields[0]]
    if actual_header != header1.split() and actual_header != header2.split():
        parser.error('First line of tension table must be "%s" or "%s" but is "%s"'
                %(header1, header2, " ".join(actual_header)))
    if not all(len(f)==len(fields[0]) for f in fields):
        parser.error('Invalid format for AFM file')
    fields = fields[1:]
    n_spring = len(fields)

    atom             = []
    spring_const     = []
    starting_tip_pos = []
    velocities       = []

    for i,f in enumerate(fields):
        res = int(f[0])
        msg = 'AFM energy specified for residue %i (zero is first residue) but there are only %i residues in the FASTA'
        if not (0 <= res < len(fasta)):
            raise ValueError(msg % (res, len(fasta)))
        atom.append( int(f[0])*3 + 1 )  # restrain the CA atom in each residue
        spring_const.append(float(f[1]))
        if len(fields[0]) == 8:
            starting_tip_pos.append( [float(x) for x in (f[2],f[3],f[4])] )
            velocities.append( [float(x) for x in (f[5],f[6],f[7])] )
        else:
            velocities.append( [float(x) for x in (f[2],f[3],f[4])] )

    atom             = np.array(atom)
    spring_const     = np.array(spring_const)
    starting_tip_pos = np.array(starting_tip_pos)
    velocities       = np.array(velocities)

    # create the moving const node
    if len(fields[0]) == 5:
        moving_node = MovingConst3DCoord2(t,  None, 'pos', time_initial, time_step, atom, velocities)
    else:
        moving_node = MovingConst3DCoord1(t,  None, 'pos', time_initial, time_step, velocities, starting_tip_pos )

    # the distance between pos and moving_node
    pairs = []
    for i in range(n_spring):
        pairs.append([atom[i], i ])
    pairs = np.array(pairs)
    pos_moving_dist  = DistanceCoord(t, 'moving', 'pos', moving_node, pairs)

    # spring potential
    dim = 0
    pbc = 0
    box_len = 0.0
    SpringPotential(t, 'pulling', pos_moving_dist, dim, pbc, box_len, np.arange(n_spring), np.zeros(n_spring), spring_const)

def write_tension(parser, fasta, tension_table):
    fields = [ln.split() for ln in open(tension_table,'r')]
    header = 'residue tension_x tension_y tension_z'
    actual_header = [x.lower() for x in fields[0]]
    if actual_header != header.split():
        parser.error('First line of tension table must be "%s" but is "%s"'
                %(header," ".join(actual_header)))
    if not all(len(f)==len(fields[0]) for f in fields):
        parser.error('Invalid format for tension file')
    fields = fields[1:]
    n_spring = len(fields)

    g = t.create_group(t.root.input.potential, 'tension')
    g._v_attrs.arguments = np.array([b'pos'])

    atom    = np.zeros((n_spring,), dtype='i')
    tension = np.zeros((n_spring,3))

    for i,f in enumerate(fields):
        res = int(f[0])
        msg = 'tension energy specified for residue %i (zero is first residue) but there are only %i residues in the FASTA'
        if not (0 <= res < len(fasta)): raise ValueError(msg % (res, len(fasta)))
        atom[i] = int(f[0])*3 + 1  # restrain the CA atom in each residue
        tension[i] = [float(x) for x in (f[1],f[2],f[3])]

    create_array(g, 'atom',    obj=atom)
    create_array(g, 'tension_coeff', obj=tension)

#=========================================================================
#                              plumed
#=========================================================================

def write_plumed(fasta, plumedFile, T=0.9, stepsize=0.009, just_print=False ):

    n_res     = len(fasta)
    n_atoms   = n_res*3

    g = t.create_group(t.root.input.potential, 'plumedforce')
    g._v_attrs.arguments = np.array([b'pos'])
    g._v_attrs.n_atoms   = n_atoms
    g._v_attrs.dt        = stepsize
    g._v_attrs.kbT       = T

    if just_print:
        g._v_attrs.just_print = 1;
    else:
        g._v_attrs.just_print = 0;

    create_array(g, 'plumedFile', obj=np.array([plumedFile]))

#=========================================================================
#                               contact 
#=========================================================================


def write_contact_energies(parser, fasta, contact_table):

    fields = [ln.split() for ln in open(contact_table,'r')]
    header_fields = 'residue1 residue2 energy distance transition_width'.split()
    if [x.lower() for x in fields[0]] != header_fields:
        parser.error('First line of contact energy table must be "%s"'%(" ".join(header_fields)))
    if not all(len(f)==len(header_fields) for f in fields):
        parser.error('Invalid format for contact file')
    fields = fields[1:]
    n_contact = len(fields)

    g = t.create_group(t.root.input.potential, 'contact')
    g._v_attrs.arguments = np.array([b'placement_fixed_point_only_CB'])

    id     = np.zeros((n_contact,2), dtype='i')
    energy = np.zeros((n_contact,))
    dist   = np.zeros((n_contact,))
    width  = np.zeros((n_contact,))

    for i,f in enumerate(fields):
        id[i] = (int(f[0]), int(f[1]))
        msg = 'Contact energy specified for residue %i (zero is first residue) but there are only %i residues in the FASTA'
        if not (0 <= id[i,0] < len(fasta)): raise ValueError(msg % (id[i,0], len(fasta)))
        if not (0 <= id[i,1] < len(fasta)): raise ValueError(msg % (id[i,1], len(fasta)))

        energy[i] = float(f[2])
        dist[i]   = float(f[3])
        width[i]  = float(f[4])  # compact_sigmoid cuts off at distance +/- width

        if width[i] <= 0.: raise ValueError('Cannot have negative contact transition_width')

    # 0-based indexing sometimes trips up users, so give them a quick check
    highlight_residues('residues that participate in any --contact potential in uppercase', fasta, id.ravel())
    if energy.max() > 0.:
        print ('\nWARNING: Some contact energies are positive (repulsive).\n'+
                 '         Please ignore this warning if you intentionally have repulsive contacts.')

    create_array(g, 'id',       obj=id)
    create_array(g, 'energy',   obj=energy)
    create_array(g, 'distance', obj=dist)
    create_array(g, 'width',    obj=width)

def write_cooperation_contacts(parser, fasta, table_list):
    table_list = table_list.split(',')
    for j, contact_table in enumerate(table_list):
        fields = [ln.split() for ln in open(contact_table,'r')]
        header_fields = 'residue1 residue2 energy distance transition_width'.split()
        if [x.lower() for x in fields[0]] != header_fields:
            parser.error('First line of contact energy table must be "%s"'%(" ".join(header_fields)))
        if not len(fields[0])==len(header_fields):
            parser.error('Invalid format for contact file')
        fields = fields[1:]
        n_contact = len(fields)

        g = t.create_group(t.root.input.potential, 'cooperation_contacts_%d'%j)
        g._v_attrs.arguments = np.array(['placement_fixed_point_only_CB'])

        id = np.zeros((n_contact,2), dtype='i')

        for i,f in enumerate(fields):
            id[i] = (int(f[0]), int(f[1]))
            msg = 'Contact energy specified for residue %i (zero is first residue) but there are only %i residues in the FASTA'
            if not (0 <= id[i,0] < len(fasta)): raise ValueError(msg % (id[i,0], len(fasta)))
            if not (0 <= id[i,1] < len(fasta)): raise ValueError(msg % (id[i,1], len(fasta)))

            if i == 0:
                energy = float(f[2])
                dist   = float(f[3])
                width  = float(f[4])  # compact_sigmoid cuts off at distance +/- width

                if width <= 0.: raise ValueError('Cannot have negative contact transition_width')

        # 0-based indexing sometimes trips up users, so give them a quick check
        highlight_residues('residues that participate in any --cooperation-contacts potential in uppercase', fasta, id.ravel())

        create_array(g, 'id',       obj=id)
        g._v_attrs.energy = energy
        g._v_attrs.dist = dist
        g._v_attrs.width = width

#=========================================================================
#                                other
#=========================================================================

def write_external_pairs_potential(pair_potential_file, use_atoms, used_percent ):
    # create a linear chain

    with tb.open_file(pair_potential_file) as epp:
        pairs              = epp.root.pairs[:]
        coefficients       = epp.root.coefficients[:]
        atoms              = epp.root.atoms[:]
        center             = epp.root.centers[:]
        atoms[atoms=='CB'] = 0
        atoms[atoms=='CA'] = 1
        atoms              = np.array(atoms, dtype=int)

    grp = t.create_group(t.root.input.potential, 'external_pairs_potential')
    grp._v_attrs.arguments = np.array([b'placement_fixed_point_only_CB', b'pos'])
   
    scales = atoms*3
    scales[atoms==0] = 1
    pairs = pairs*scales+atoms

    if use_atoms == "CA":
        ndx = np.where((atoms[:,0] == 1) * (atoms[:,1] == 1))
    elif use_atoms == "CB":
        ndx = np.where((atoms[:,0] == 0) + (atoms[:,1] == 0))
    elif use_atoms == "ALL":
        ndx = np.where((atoms[:,0] != 2) + (atoms[:,1] != 2))

    pair = pairs[ndx]
    atom = atoms[ndx]
    coefficient = coefficients[ndx]
        
    create_array(grp, 'pair',        obj=pair)
    create_array(grp, 'atom',        obj=atom)
    create_array(grp, 'coefficient', obj=coefficient)
    create_array(grp, 'center',      obj=center)
    create_array(grp, 'used',        obj=np.array([used_percent]))

def write_sidechain_radial(fasta, library, excluded_residues, suffix=''):
    g = t.create_group(t.root.input.potential, 'radial'+suffix)
    g._v_attrs.arguments = np.array([b'placement_fixed_point_only_CB'])
    for res_num in excluded_residues:
        if not (0<=res_num<len(fasta)):
            raise ValueError('Residue number %i is invalid'%res_num)

    residues = sorted(set(np.arange(len(fasta))).difference(excluded_residues))

    with tb.open_file(library) as params:
        resname2restype = dict((x,i) for i,x in enumerate(params.root.names[:]))
        n_type = len(resname2restype)

        create_array(g, 'index', obj=np.array(residues))
        create_array(g, 'type',  obj=np.array([resname2restype[x] for x in fasta[residues]]))
        create_array(g, 'id',    obj=np.array(residues))  # FIXME update for chain breaks
        create_array(g, 'interaction_param', obj=params.root.interaction_param[:])


def parse_segments(s):
    ''' Parse segments of the form 10-30,50-60 '''
    import argparse
    import re

    if re.match('^([0-9]+(-[0-9]+)?)(,[0-9]+(-[0-9]+)?)*$', s) is None:
        raise argparse.ArgumentTypeError('segments must be of the form 10-30,45,72-76 or similar')

    def parse_seg(x):
        atoms = x.split('-')
        if len(atoms) == 1:
            return np.array([int(atoms[0])])
        elif len(atoms) == 2:
            return np.arange(int(atoms[0]),1+int(atoms[1]))  # inclusive on both ends
        else:
            raise RuntimeError('the impossible happened.  oops.')

    ints = np.concatenate([parse_seg(a) for a in s.split(',')])
    ints = np.array(sorted(set(ints)))   # remove duplicates and sort
    return ints


def parse_float_pair(s):
    import argparse
    import re

    args = s.split(',')
    if len(args) != 2:
        raise argparse.ArgumentTypeError('must be in the form -2.0,-1.0 or similar (exactly 2 numbers)')

    return (float(args[0]), float(args[1]))


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Prepare input file',
            usage='use "%(prog)s --help" for more information')

    parser.add_argument('--config', required=True,
            help='[required] Upside config file')

    parser.add_argument('--group', default='',
            help='Table of group.  Each line must contain 1 fields and the first line '+
            'must be "residues". The following lines is the residue range like this: "2-10,12".')

    parser_grp1 = parser.add_mutually_exclusive_group()
    parser_grp1.add_argument('--cavity-radius', default=0., type=float,
            help='Enclose the whole simulation in a radial cavity centered at the origin to achieve finite concentration '+
            'of protein.  Necessary for multichain simulation.')
    parser_grp1.add_argument('--heuristic-cavity-radius', default=0., type=float,
        help='Set the cavity radius to this provided scale factor times the max distance between com\'s and atoms of the chains.')
    parser_grp1.add_argument('--cavity-radius-from-config', default='', help='Config file with cavity radius set. Useful for applying'+
            ' the same heuristic cavity of bound complex config to unbound counterpart')
    parser.add_argument('--make-unbound', action='store_true',
        help='Separate chains into different corners of a cavity that you set with one of the cavity options.')

    parser.add_argument('--fixed-wall', default='',
            help='Table of fixed wall springs.  Each line must contain 5-7 fields and the first line '+
            'must contain one of the following headers: ' +
            '"residue radius spring_const wall_type x0", ' +
            '"residue radius spring_const wall_type y0", ' +
            '"residue radius spring_const wall_type z0", ' +
            '"residue radius spring_const wall_type x0 y0", ' +
            '"residue radius spring_const wall_type x0 z0", ' +
            '"residue radius spring_const wall_type y0 z0", ' +
            '"residue radius spring_const wall_type x0 y0 z0". ' + 
            'Depending on the header, 1-D, 2-D or 3-D wall spring restraint force is added when the distance between the CA atom of given residue ' + 
            '(starting from 0) and the reference coordinate is lower (wall_type=0) or higher (wall_type=1) than a preset length. ' +
            'The "residue" column supports you to specify the residue range like this: "2-10,12". These residues will share the same reference ' +
            'position and the wall spring properties')

    parser.add_argument('--pair-wall', default='',
            help='Table of wall springs for the distance between two atoms.  The first line '+
            'must contain one of the following headers: ' +
            '"residue1 residue2 radius spring_const wall_type x", ' +
            '"residue1 residue2 radius spring_const wall_type y", ' +
            '"residue1 residue2 radius spring_const wall_type z", ' +
            '"residue1 residue2 radius spring_const wall_type xy", ' +
            '"residue1 residue2 radius spring_const wall_type xz", ' +
            '"residue1 residue2 radius spring_const wall_type yz", ' +
            '"residue1 residue2 radius spring_const wall_type xyz". ' + 
            'The last column is used to specify the dimension of the wall spring potential.' + 
            'The second and subsequent lines must contain 5 fields to define the residue pair, upper limit (or lower limit), ' + 
            'spring constant and wall type (0: lower wall, 1: upper wall) respectively.' +
            'The wall works when the distance between two CA atoms of given residues ' + 
            'is lower (wall_type=0) or higher (wall_type=1) than a preset length.')

    parser.add_argument('--fixed-spring', default='',
            help='Table of fixed springs.  Each line must contain 4-6 fields and the first line '+
            'must contain one of the following headers: ' +
            '"residue radius spring_const x0", ' +
            '"residue radius spring_const y0", ' +
            '"residue radius spring_const z0", ' +
            '"residue radius spring_const x0 y0", ' +
            '"residue radius spring_const x0 z0", ' +
            '"residue radius spring_const y0 z0", ' +
            '"residue radius spring_const x0 y0 z0". ' + 
            'Depending on the header, 1-D, 2-D or 3-D spring restraint force is added when the distance between the CA atom of given residue ' + 
            '(starting from 0) and the reference coordinate is different to a preset length "radius". ' +
            'The "residue" column supports you to specify the residue range like this: "2-10,12". These residues will share the same reference ' +
            'position and the spring properties')

    parser.add_argument('--group-const-spring', default='',
            help='Table of fixed springs.  Each line must contain 4-6 fields and the first line '+
            'must contain one of the following headers: ' +
            '"group radius spring_const x0", ' +
            '"group radius spring_const y0", ' +
            '"group radius spring_const z0", ' +
            '"group radius spring_const x0 y0", ' +
            '"group radius spring_const x0 z0", ' +
            '"group radius spring_const y0 z0", ' +
            '"group radius spring_const x0 y0 z0". ' + 
            'Depending on the header, 1-D, 2-D or 3-D spring restraint force is added when the distance between the group center ' + 
            'and the reference coordinate is different to a preset length "radius". ' )


    parser.add_argument('--pair-spring', default='',
            help='Table of springs for the distance between two atoms.  The first line '+
            'must contain one of the following headers: ' +
            '"residue1 residue2 radius spring_const x", ' +
            '"residue1 residue2 radius spring_const y", ' +
            '"residue1 residue2 radius spring_const z", ' +
            '"residue1 residue2 radius spring_const xy", ' +
            '"residue1 residue2 radius spring_const xz", ' +
            '"residue1 residue2 radius spring_const yz", ' +
            '"residue1 residue2 radius spring_const xyz". ' + 
            'The last column is used to specify the dimension of the spring potential.' + 
            'The second and subsequent lines must contain 4 fields to define the residue pair, the balance position, ' + 
            'and the spring constant respectively.')

    parser.add_argument('--group-dist-spring', default='',
            help='Table of springs for the distance between two atoms.  The first line '+
            'must contain one of the following headers: ' +
            '"group1 group2 radius spring_const x", ' +
            '"group1 group2 radius spring_const y", ' +
            '"group1 group2 radius spring_const z", ' +
            '"group1 group2 radius spring_const xy", ' +
            '"group1 group2 radius spring_const xz", ' +
            '"group1 group2 radius spring_const yz", ' +
            '"group1 group2 radius spring_const xyz". ' + 
            'The last column is used to specify the dimension of the spring potential.' + 
            'The second and subsequent lines must contain 4 fields to define the residue pair, the balance position, ' + 
            'and the spring constant respectively.')
    
    parser.add_argument('--nail-spring', default='',
            help='Table of nail springs to nail the CA atoms to the initial position. The first line '+
            'must contain one of the following headers: ' +
            '"residue spring_const x", ' +
            '"residue spring_const y", ' +
            '"residue spring_const z", ' +
            '"residue spring_const xy", ' +
            '"residue spring_const xz", ' +
            '"residue spring_const yz", ' +
            '"residue spring_const xyz". ' + 
            'The last element of the header is used to specify the dimension of the spring potential.' + 
            'The second and subsequent lines must contain 2 fields to define the residue and ' + 
            'the spring constant respectively. You to specify the residue range like this: "2-10,12".')

    parser.add_argument('--group-nail-spring', default='',
            help='Table of nail springs to nail the CA atoms to the initial position. The first line '+
            'must contain one of the following headers: ' +
            '"group spring_const x", ' +
            '"group spring_const y", ' +
            '"group spring_const z", ' +
            '"group spring_const xy", ' +
            '"group spring_const xz", ' +
            '"group spring_const yz", ' +
            '"group spring_const xyz". ' + 
            'The last element of the header is used to specify the dimension of the spring potential.' + 
            'The second and subsequent lines must contain 2 fields to define the group id and ' + 
            'the spring constant respectively.')

    parser.add_argument('--restraint-group', default=[], action='append', type=parse_segments,
            help='List of residues in the protein.  The residue list should be of a form like ' +
            '--restraint-group 10-13,17,19-21 and that list would specify all the atoms in '+
            'residues 10,11,12,13,17,19,20,21. '+
            'Each atom in the specified residues will be randomly connected to atoms in other residues by ' +
            'springs with equilibrium distance given by the distance of the atoms in the initial structure.  ' +
            'Multiple restraint groups may be specified by giving the --restraint-group flag multiple times '
            'with different residue lists.  The strength of the restraint is given by --restraint-spring-constant')

    parser.add_argument('--apply-restraint-group-to-each-chain', action='store_true',
            help='Use indices of chain first residues recorded during PDB_to_initial_structure to automate'+
            ' --restraint-group for chains. Requires --chain-break-from-file.')

    parser.add_argument('--restraint-spring-constant', default=4., type=float,
            help='Spring constant used to restrain atoms in a restraint group (default 4.) ')

    parser.add_argument('--ask-before-using-AFM', default='',
            help='Table of tip positions and pulling velocitis for mimicing AFM pulling experiment in the constant velocity mode. ' +
            'Each line must contain 8 fields and the first line must contain ' +
            '"residue spring_const tip_pos_x tip_pos_y tip_pos_z pulling_vel_x pulling_vel_y pulling_vel_z". ' +
            'The residue will be pulled in the direction (pulling_vel_x, pulling_vel_y, pulling_vel_z) by its CA atom, ' +
            'which is attached to the tip at (tip_pos_x, tip_pos_y, tip_pos_z). ' +
            'The magnitude of the pulling velocity vector sets the pulling speed. The unit is: angstrom/time_step. ' +
            'The spring_const is in the unit of kT/angstrom^2. At T = 298.15 K, it equals 41.14 pN/angstrom. ' + 
            'Note: consult with the developer before using this AFM function. '
            'Notice that to use this option, you have to either disable '+
            'the recentering (soluble protein) or distable z recentering for membrane protein. ' +
            'To do that, you will need to add --disable-recentering or --disable-z-recentering in UPSIDE EXECUTABLE arguments '+
            '(not arguments for upside_config.py or advanced_config.py)')
    parser.add_argument('--AFM-time-initial', default=0., type=float,
            help='Time initial for AFM pulling simulation. The default value is 0. ' +
            'WARNING: do not change this value unless the simulation is a continuation of a previous one. ' +
            'To set the time initial, check the /root/output/time_estimate in the output h5 file. ' )
    parser.add_argument('--AFM-time-step', default=0.009, type=float,
            help='Time step for AFM pulling simulation. The default value is 0.009. ' +
            'WARNING: this should be the same as the global time step, which is set to 0.009 by default. Change this value accordingly.')

    parser.add_argument('--tension', default='',
            help='Table of linear tensions.  Each line must contain 4 fields and the first line '+
            'must contain "residue tension_x tension_y tension_z".  The residue will be pulled in the '+
            'direction (tension_x,tension_y,tension_z) by its CA atom.  The magnitude of the tension vector '+
            'sets the force.  Units are kT/Angstrom. Notice that to use this option, you have to either disable '+
            'the recentering (soluble protein) or distable z recentering for membrane protein. ' +
            'To do that, you will need to add --disable-recentering or --disable-z-recentering in UPSIDE EXECUTABLE arguments '+
            '(not arguments for upside_config.py or advanced_config.py)')

    parser.add_argument('--contact-energies', default='',
            help='Path to text file that defines a contact energy function.  The first line of the file should ' +
            'be a header containing "residue1 residue2 energy distance transition_width", and the remaining '+
            'lines should contain space separated values.  The form of the interaction is approximately '+
            'sigmoidal but the potential is constant outside (distance-transition_width,distance+transition_width).'+
            '  This potential is approximately twice as sharp as a standard sigmoid with the same width as the '+
            'specified transition_width.  The location x_residue is approximately the CB position of the '+
            'residue.')
    parser.add_argument('--cooperation-contacts', default='',
            help='Path to text file that defines a cooperation contacts function. The first line of the file ' +
            'should be a header containing "residue1 residue2 energy distance transition_width", and the remaining ' +
            'lines should contain space separated values. Only the values for energy, distance, and ' +
            'transition_width from the first entry are used. It is calculated as a product of approximate unit ' +
            'sigmoidal functions, scaled by the energy param. See help for --contact-energies for details on the ' +
            'form of the individual sigmoidal functions')

    parser.add_argument('--plumed', default='', help='Path to text file that serves as a standard control file to configure plumed simulation.')
    parser.add_argument('--plumed-time-step', default=0.009, type=float,
            help='set time step to PLUMED. The default value is 0.009. ' +
            'WARNING: this should be the same as the global time step, which is set to 0.009 by default. Change this value accordingly.')
    parser.add_argument('--plumed-temperature', default=0.9, type=float,
            help='set temperature to PLUMED. The default value is 0.9. ' +
            'WARNING: this should be the same as the global temperature, which is set by --temperature.')

    parser.add_argument('--external-pairs-table-potential', default='', help='User-defined table-type pair interactions. (pair and potential)')
    parser.add_argument('--external-pairs-type',            default='', help='User-defined table-type pair interactions. (used atoms: CA, CB or ALL)')
    parser.add_argument('--external-pairs-used-percent',    default=1.0, type=float, help='User-defined table-type pair interactions. (used percent)')

    args = parser.parse_args()


    global n_atom, t, potential
    t = tb.open_file(args.config,'a')
    potential = t.root.input.potential
    fasta     = t.root.input.sequence[:]
    pos       = t.root.input.pos[:]
    n_res     = len(fasta)
    n_atom    = 3*n_res

    use_append_const3D_to_pos = False

    multi_chain = False
    has_rl_info = False
    require_backbone_point = False
    require_affine = False

    # Record of adv config options
    args_group = t.create_group(t.root.input, 'adv_args')
    for k,v in sorted(vars(args).items()):
        # h5 has header message size limit, _v_attrs values fall under this.
        # restraint_group values for large proteins can exceed this as fully expanded array, so
        # turn them back into compact range form, i.e. ["0-3"], instead of [[0, 1, 2, 3]] 
        if k == "restraint_group":
            v_compact = []
            for restr_grp in v:
                G=(list(x) for _,x in groupby(restr_grp, lambda x,c=count(): next(c)-x))
                v_compact.append(",".join("-".join(map(str,(g[0],g[-1])[:len(g)])) for g in G))
            v = v_compact
        args_group._v_attrs[k] = v
    args_group._v_attrs['invocation'] = ' '.join(sys.argv[:])


    if 'chain_break' in t.root.input:
        chain_first_residue = t.root.input.chain_break.chain_first_residue[:]
        multi_chain = True
        n_chains = len(chain_first_residue) + 1
        if 'rl_chains' in t.root.input.chain_break:
            has_rl_info = True
            rl_chains = t.root.input.chain_break.rl_chains[:]

    if args.apply_restraint_group_to_each_chain and ('chain_break' not in t.root.input): 
        parser.error('--apply-restraint-group-to-each-chain in advanced_config requires --chain-break-from-file in upside_config')

    if args.heuristic_cavity_radius:
        if n_chains < 2:
            eprint('WARNING: --heuristic-cavity-radius requires at least 2 chains. Skipping setting up cavity')
        else:
            if not has_rl_info: # all chains separately
                com_list = []
                for i in range(n_chains):
                    first_res, next_first_res = chain_endpts(n_res, chain_first_residue, i)
                    com_list.append(pos[first_res*3:next_first_res*3,:,0].mean(axis=0))
            else: # receptor and ligand chains considered as groups
                # receptor com
                first_res = chain_endpts(n_res, chain_first_residue, 0)[0]
                next_first_res = chain_endpts(n_res, chain_first_residue, rl_chains[0]-1)[1]
                r_com = pos[first_res*3:next_first_res*3,:,0].mean(axis=0)

                # ligand com
                first_res = chain_endpts(n_res, chain_first_residue, rl_chains[0])[0]
                next_first_res = chain_endpts(n_res, chain_first_residue, n_chains-1)[1]
                l_com = pos[first_res*3:next_first_res*3,:,0].mean(axis=0)

                com_list = [r_com, l_com]
            # Distance between chain com and all atoms
            com_dist_list = []
            for i in range(len(com_list)):
                for j in range(n_atom):
                        com_dist_list.append(vmag(com_list[i]-pos[j,:,0]) + 2.)

            args.cavity_radius = args.heuristic_cavity_radius*max(com_dist_list)
            print ()
            print ("cavity_radius")
            print (args.cavity_radius)

    if args.cavity_radius_from_config:
        if n_chains < 2:
            eprint('WARNING: --cavity-radius-from-config requires at least 2 chains. Skipping setting up cavity')
        elif args.heuristic_cavity_radius:
            eprint('WARNING: Overwriting heuristic cavity with the one from --cavity-radius-from-config')
        else:
            t_cavity = tb.open_file(args.cavity_radius_from_config,'r')
            args.cavity_radius = t_cavity.root.input.potential.cavity_radial.radius[0]
            t_cavity.close()
            print ()
            print ("cavity_radius")
            print (args.cavity_radius)

    if args.make_unbound:
        if n_chains < 2 or n_chains > 8:
            eprint('WARNING: --make-unbound requires at least 2 and no more than 8 chains. Skipping separating chains.\n'
                   'See --record-chain-breaks in PDB_to_initial_structure.py and --chain-break-from-file in upside_config.py')
        elif not args.cavity_radius:
            eprint('WARNING: --make-unbound requires setting a cavity radius. Skipping separating chains')
        else:
            print ()
            print ("making unbound")
            pos = t.root.input.pos[:]
            displacement = np.array([[-1.,0.,0.], [1.,0.,0.],
                                     [0.,-1.,0.], [0.,1.,0.],
                                     [0.,0.,-1.], [0.,0.,1.],])
            if not has_rl_info: # separate all chains
                for j in range(n_chains):
                    first_res, next_first_res = chain_endpts(n_res, chain_first_residue, j)
                    #com = pos[first_res*3:next_first_res*3,:,0].mean(axis=0)
                    pos[first_res*3:next_first_res*3,:,0] = (pos[first_res*3:next_first_res*3,:,0] +
                            displacement[j]*0.5*args.cavity_radius) #- displacement[j]*com
            else: # keep receptor and ligand chains together
                # move receptor chains
                first_res = chain_endpts(n_res, chain_first_residue, 0)[0]
                next_first_res = chain_endpts(n_res, chain_first_residue, rl_chains[0]-1)[1]
                pick_disp = np.random.choice([0, 2, 4])
                pos[first_res*3:next_first_res*3,:,0] = pos[first_res*3:next_first_res*3,:,0] + displacement[pick_disp]*0.5*args.cavity_radius

                # move ligand chains
                first_res = chain_endpts(n_res, chain_first_residue, rl_chains[0])[0]
                next_first_res = chain_endpts(n_res, chain_first_residue, n_chains-1)[1]
                pick_disp = np.random.choice([1, 3, 5])
                pos[first_res*3:next_first_res*3,:,0] = pos[first_res*3:next_first_res*3,:,0] + displacement[pick_disp]*0.5*args.cavity_radius
            t.root.input.pos[:] = pos

    if args.group:
        write_group_node(parser, fasta, args.group)

    if args.cavity_radius:
        write_cavity_radial(args.cavity_radius)

    if args.fixed_wall:
        write_const_wall_spring(parser, fasta, args.fixed_wall)

    if args.pair_wall:
        write_pair_wall_spring(parser, fasta, args.pair_wall)

    if args.fixed_spring:
        write_const_spring(parser, fasta, args.fixed_spring)

    if args.pair_spring:
        write_pair_spring(parser, fasta, args.pair_spring)

    if args.nail_spring:
        write_nail_spring(parser, fasta, args.nail_spring)

    if args.group_dist_spring:
        write_group_dist_spring(parser, fasta, args.group_dist_spring)

    if args.group_nail_spring:
        write_group_nail_spring(parser, fasta, args.group_nail_spring)

    if args.group_const_spring:
        write_group_const_spring(parser, fasta, args.group_const_spring)

    if args.apply_restraint_group_to_each_chain and n_chains > 1:
        if has_rl_info:
            # receptor chains
            first_res      = chain_endpts(n_res, chain_first_residue, 0)[0]
            next_first_res = chain_endpts(n_res, chain_first_residue, rl_chains[0]-1)[1]
            args.restraint_group.append(np.arange(first_res, next_first_res))

            # ligand chains
            first_res = chain_endpts(n_res, chain_first_residue, rl_chains[0])[0]
            next_first_res = chain_endpts(n_res, chain_first_residue, n_chains-1)[1]
        print ("restraint_group")
        print (args.restraint_group)

    if args.restraint_group:
        print ()
        print ('Restraint groups (uppercase letters are restrained residues)')
        fasta_one_letter = ''.join(one_letter_aa[x.decode('ASCII')] for x in fasta)
        print ()
        print ("Restraint spring constant: {}".format(args.restraint_spring_constant))

        for i,restrained_residues in enumerate(args.restraint_group):
            assert np.amax(list(restrained_residues)) < len(fasta)
            highlight_residues('group_%i'%i, fasta, restrained_residues)
            print (i, restrained_residues)
            make_restraint_group(i, set(restrained_residues), pos[:,:,0], args.restraint_spring_constant)

    if args.tension and args.ask_before_using_AFM:
        print ('Nope, you cannot pull the protein using two modes. Choose one.')
    elif args.tension and not args.ask_before_using_AFM:
        write_tension(parser, fasta, args.tension)
    elif args.ask_before_using_AFM and not args.tension:
        write_pulling(parser, fasta, args.ask_before_using_AFM, args.AFM_time_initial, args.AFM_time_step)

    if args.plumed:
        write_plumed(fasta, args.plumed, args.plumed_temperature, args.plumed_time_step)

    if args.contact_energies:
        require_backbone_point = True
        write_contact_energies(parser, fasta, args.contact_energies)

    if args.cooperation_contacts:
        require_backbone_point = True
        write_cooperation_contacts(parser, fasta, args.cooperation_contacts)

    if args.external_pairs_table_potential:
        require_backbone_point = True
        if args.external_pairs_type is None:
            parser.error('--external-pairs-table-potential requires --external-pairs-type')
        else:
            atom_type = args.external_pairs_type
            if atom_type != 'CA' and atom_type != 'CB' and atom_type != 'ALL':
                parser.error('--external-pairs-type must be CA, CB or ALL')
            write_external_pairs_potential(args.external_pairs_table_potential, atom_type, args.external_pairs_used_percent)

    if use_append_const3D_to_pos:
        append_const3D_to_pos()

    if require_backbone_point:
        require_affine = True
        if not 'placement_fixed_point_only_CB' in t.root.input.potential:
            write_CB(fasta)

    if require_affine:
        if not 'affine_alignment' in t.root.input.potential:
            write_affine_alignment(n_res)

    t.close()

if __name__ == '__main__':
    main()
