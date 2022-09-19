#!/usr/bin/env python
import numpy as np
import tables as tb
import sys,os
import _pickle as cPickle
from upside_nodes import *

#---------------------------------------------------------------------------
#                     Predefined variables, types
#---------------------------------------------------------------------------

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def bstring(string):
    return bytes(string, encoding="ascii")

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

backbone_group0 = { 'ALA': 0,  'ARG': 1,  'ASN': 2,  'ASP': 3,
                    'CYS': 4,  'GLN': 5,  'GLU': 6,  'GLY': 7,
                    'HIS': 8,  'ILE': 9,  'LEU': 10, 'LYS': 11,
                    'MET': 12, 'PHE': 13, 'PRO': 19, 'SER': 14,
                    'THR': 15, 'TRP': 16, 'TYR': 17, 'VAL': 18  }

backbone_group1 = { 'ALA': 0,  'ARG': 0,  'ASN': 1,  'ASP': 1,
                    'CYS': 0,  'GLN': 0,  'GLU': 0,  'GLY': 2,
                    'HIS': 0,  'ILE': 3,  'LEU': 0,  'LYS': 0,
                    'MET': 0,  'PHE': 0,  'PRO': 5,  'SER': 0,
                    'THR': 4,  'TRP': 0,  'TYR': 0,  'VAL': 3  }

backbone_group2 = { 'ALA': 0,  'ARG': 0,  'ASN': 1,  'ASP': 1,
                    'CYS': 2,  'GLN': 0,  'GLU': 0,  'GLY': 3,
                    'HIS': 0,  'ILE': 4,  'LEU': 0,  'LYS': 0,
                    'MET': 2,  'PHE': 0,  'PRO': 7,  'SER': 5,
                    'THR': 6,  'TRP': 2,  'TYR': 0,  'VAL': 4  }

backbone_group4 = { 'ALA': 0,  'ARG': 0,  'ASN': 0,  'ASP': 0,
                    'CYS': 0,  'GLN': 0,  'GLU': 0,  'GLY': 1,
                    'HIS': 0,  'ILE': 0,  'LEU': 0,  'LYS': 0,
                    'MET': 0,  'PHE': 0,  'PRO': 2,  'SER': 0,
                    'THR': 0,  'TRP': 0,  'TYR': 0,  'VAL': 0  }

backbone_group5 = { 'ALA': 0,  'ARG': 0,  'ASN': 0,  'ASP': 0,
                    'CYS': 0,  'GLN': 0,  'GLU': 0,  'GLY': 0,
                    'HIS': 0,  'ILE': 0,  'LEU': 0,  'LYS': 0,
                    'MET': 0,  'PHE': 0,  'PRO': 0,  'SER': 0,
                    'THR': 0,  'TRP': 0,  'TYR': 0,  'VAL': 0  }

#---------------------------------------------------------------------------
#                            Utility functions
#---------------------------------------------------------------------------

def highlight_residues(name, fasta, residues_to_highlight):
    fasta_one_letter = [one_letter_aa[x] for x in fasta]
    residues_to_highlight = set(residues_to_highlight)
    print ('%s:  %s' % (name, ''.join((f.upper() if i in residues_to_highlight else f.lower()) for i,f in enumerate(fasta_one_letter))))

def vmag(x):
    assert x.shape[-1] == 3
    return np.sqrt(x[...,0]**2+x[...,1]**2+x[...,2]**2)

def create_array(grp, nm, obj=None):
    return t.create_earray(grp, nm, obj=obj, filters=default_filter)


#-------------------------------------
# build the random initial structure

def make_tab_matrices(phi, theta, bond_length):
    '''TAB matrices are torsion-angle-bond affine transformation matrices'''
    phi         = np.asarray(phi)
    theta       = np.asarray(theta)
    bond_length = np.asarray(bond_length)

    assert phi.shape == theta.shape == bond_length.shape
    r = np.zeros(phi.shape + (4,4), dtype=(phi+theta+bond_length).dtype)

    cp = np.cos(phi  ); sp = np.sin(phi  )
    ct = np.cos(theta); st = np.sin(theta)
    l  = bond_length

    r[...,0,0]=   -ct; r[...,0,1]=    -st; r[...,0,2]=   0; r[...,0,3]=   -l*ct;
    r[...,1,0]= cp*st; r[...,1,1]= -cp*ct; r[...,1,2]= -sp; r[...,1,3]= l*cp*st;
    r[...,2,0]= sp*st; r[...,2,1]= -sp*ct; r[...,2,2]=  cp; r[...,2,3]= l*sp*st;
    r[...,3,0]=     0; r[...,3,1]=      0; r[...,3,2]=   0; r[...,3,3]=       1;

    return r

def construct_equilibrium_structure(rama, angles, bond_lengths):
    assert rama.shape == angles.shape == bond_lengths.shape
    n_res = rama.shape[0]
    n_atom = 3*n_res
    assert rama.shape == (n_res,3)

    t = np.zeros(n_atom)
    a = angles.ravel()
    b = bond_lengths.ravel()

    t[3::3] = rama[:-1,1]
    t[4::3] = rama[:-1,2]
    t[5::3] = rama[1: ,0]

    transforms = make_tab_matrices(t,a,b)
    curr_affine = np.eye(4)
    pos = np.zeros((3*n_res,3))

    # right apply all transformations

    for i,mat in enumerate(transforms):
        curr_affine = np.dot(curr_affine, mat)
        pos[i] = curr_affine[:3,3]
    return pos


def random_initial_config(n_res):
    # a reasonable model where the chain grows obeying sensible angles and omegas
    rama    = np.random.random((n_res,3))*2*np.pi - np.pi
    angles  = np.zeros_like(rama)
    lengths = np.zeros_like(rama)

    rama[:,2] = np.pi   # all trans omega's

    angles[:,0] = 120.0*deg  # CA->C->N angle
    angles[:,1] = 120.0*deg  # C->N->CA angle
    angles[:,2] = 109.5*deg  # N->CA->C angle

    lengths[:,0] = 1.453
    lengths[:,1] = 1.526
    lengths[:,2] = 1.300
    return construct_equilibrium_structure(rama, angles, lengths)
#-------------------------------------

def compact_sigmoid(x, sharpness):
    y = x*sharpness;
    result = 0.25 * (y+2) * (y-1)**2
    result = np.where((y< 1), result, np.zeros_like(result))
    result = np.where((y>-1), result, np.ones_like (result))
    return result

def double_compact_sigmoid(x, half_width, sharpness):
    return compact_sigmoid(x-half_width, sharpness) * compact_sigmoid(-x-half_width, sharpness)

def angular_compact_double_sigmoid(theta, center, half_width, sharpness):
    dev = theta-center
    dev = np.where((dev< np.pi), dev, dev-2*np.pi)
    dev = np.where((dev>-np.pi), dev, dev+2*np.pi)
    return double_compact_sigmoid(dev, half_width, sharpness)

def rama_box(rama, center, half_width, sharpness):

    assert rama.shape[-1] == center.shape[-1] == half_width.shape[-1] == 2

    s = center.shape[:-1]
    if not s:
        return (angular_compact_double_sigmoid(rama[...,0], rama, center, half_width, sharpness)*
                angular_compact_double_sigmoid(rama[...,1], rama, center, half_width, sharpness))
    else:
        result = np.zeros(rama.shape[:-1] + center.shape[:-1])
        for inds in np.indices(s).reshape((len(s),-1)).T:
            inds = tuple(inds)
            if len(inds) == 1: inds = inds[0]
            value = (
                    angular_compact_double_sigmoid(rama[...,0], center[inds,0], half_width[inds,0], sharpness)*
                    angular_compact_double_sigmoid(rama[...,1], center[inds,1], half_width[inds,1], sharpness))
            result[...,inds] = value
        return result

def read_fasta(file_obj):
    lines = list(file_obj)
    assert lines[0][0] == '>'
    one_letter_seq = ''.join(x.strip().replace('\r','') for x in lines[1:])
    seq = []
    cis_state = False
    for a in one_letter_seq:
        if cis_state:
            assert a == 'P'  # proline must follow start
            seq.append('CPR')
            cis_state = False
        elif a == "*":
            cis_state = True
        else:
            seq.append(three_letter_aa[a])
    return np.array(seq)

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

def find_multichain_terms(ids, chain_starts):
    assert len(ids.shape) == 2
    chain_starts = np.array(chain_starts, dtype='i')
    assert len(chain_starts.shape) == 1
    chain_num = (ids[:][:,:,None] >= chain_starts[None,None,:]).sum(axis=-1)
    multichain = np.array([len(set(x)) <= 1 for x in chain_num])
    return multichain

#---------------------------------------------------------------------------
#                                CoordNodes
#---------------------------------------------------------------------------

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

def write_infer_H_O(fasta, excluded_residues):
    n_res = len(fasta)
    # note that proline is not an hbond donor since it has no NH
    donor_residues    = np.array([i for i in range(n_res) if i>0       and i not in excluded_residues and fasta[i]!='PRO'])
    acceptor_residues = np.array([i for i in range(n_res) if i<n_res-1 and i not in excluded_residues])

    print()
    print ('hbond, %i donors, %i acceptors in sequence' % (len(donor_residues), len(acceptor_residues)))

    H_bond_length = 0.88
    O_bond_length = 1.24

    grp = t.create_group(potential, 'infer_H_O')
    grp._v_attrs.arguments = np.array([b'pos'])

    donors    = t.create_group(grp, 'donors')
    acceptors = t.create_group(grp, 'acceptors')

    create_array(donors,    'residue', obj=donor_residues)
    create_array(acceptors, 'residue', obj=acceptor_residues)

    create_array(donors,    'bond_length', obj=H_bond_length*np.ones(len(   donor_residues)))
    create_array(acceptors, 'bond_length', obj=O_bond_length*np.ones(len(acceptor_residues)))

    d_id = np.array((-1,0,1))[None,:] + 3*donor_residues   [:,None]
    a_id = np.array(( 1,2,3))[None,:] + 3*acceptor_residues[:,None]
    create_array(donors,    'id', obj=d_id)
    create_array(acceptors, 'id', obj=a_id)

    if n_chains > 1:
        print ('Checking %s' % "infer_H_O")
        print ()

        is_bad =(any((1-find_multichain_terms(d_id, chain_starts))) or any((1-find_multichain_terms(a_id, chain_starts))))
        if is_bad:
            print ('!!! Error: You must use --hbond-excluded-residues in upside_config to produce chain breaks for the hbonds !!!')
            print ()

def write_weighted_pos(sc_node_name, pl_node_name): 
    # Bring position and probability together for the side chains
    wgrp = t.create_group(potential, 'weighted_pos')
    wgrp._v_attrs.arguments = np.array([bstring(sc_node_name), bstring(pl_node_name)])
    sc_node = t.get_node(t.root.input.potential, sc_node_name)
    n_sc = sc_node.affine_residue.shape[0]
    create_array(wgrp, 'index_pos',   np.arange(n_sc))
    create_array(wgrp, 'index_weight', np.arange(n_sc))

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

def write_rama_coord():
    grp = t.create_group(potential, 'rama_coord')
    grp._v_attrs.arguments = np.array([b'pos'])
    n_res = n_atom/3
    N_id = 3*np.arange(n_res)
    id = np.column_stack((N_id-1,N_id,N_id+1,N_id+2,N_id+3))
    id[id>=n_atom] = -1  # last atom is non-existent (as is first)
                         #   and non-existence is indicated by -1
    create_array(grp, 'id', id)

def write_rama_coord2():
    n_res = n_atom/3
    N_id = 3*np.arange(n_res)
    id1 = np.column_stack((N_id-1,N_id,N_id+1,N_id+2))
    id2 = np.column_stack((N_id,N_id+1,N_id+2,N_id+3))
    id2[id2>=n_atom] = -1  # last atom is non-existent (as is first)
    id3 = np.column_stack((N_id-1,N_id,N_id+1,N_id+2,N_id+3))
    id3[id3>=n_atom] = -1  # last atom is non-existent (as is first)

    if n_chains > 1:
        multichain_locs = (1-find_multichain_terms(id3, chain_starts)).nonzero()[0]
        print ()
        print ('Editing %s (%i rows)' %('rama_coord', len(multichain_locs)))
        for loc in multichain_locs:
            ids = id3[loc]
            assert ids.shape == (5,)
            chain_num = (ids[:,None]>=chain_starts).sum(axis=-1)
            if not (chain_num[1] == chain_num[2] == chain_num[3] and (chain_num[0]==chain_num[1] or chain_num[3]==chain_num[4])): 
                raise ValueError("Weird rama_coord %i, unable to proceed" % loc)

            if chain_num[0]==chain_num[1]: 
                id3[loc,4] = -1  # cut psi
                id2[loc,3] = -1  # cut psi
            else: 
                id3[loc,0] = -1  # cut phi
                id1[loc,0] = -1  # cut phi

    phi_grp = t.create_group(potential, 'Dihedral_phi')
    phi_grp._v_attrs.arguments = np.array([b'pos'])
    create_array(phi_grp, 'id', obj=id1)

    psi_grp = t.create_group(potential, 'Dihedral_psi')
    psi_grp._v_attrs.arguments = np.array([b'pos'])
    create_array(psi_grp, 'id', obj=id2)

    grp = t.create_group(potential, 'rama_coord')
    grp._v_attrs.arguments = np.array([b'Dihedral_phi', b'Dihedral_psi'])
    create_array(grp, 'index_pos1', np.arange(id1.shape[0]))
    create_array(grp, 'index_pos2', np.arange(id2.shape[0]))
    create_array(grp, 'id', id3)

#---------------------------------------------------------------------------
#                               Bonded terms
#---------------------------------------------------------------------------

# write dist_spring potential
def write_dist_spring(args):
    # create a linear chain
    dgrp = t.create_group(potential, 'Distance3D')
    dgrp._v_attrs.arguments = np.array([b'pos', b'pos'])

    id = np.arange(n_atom-1)
    id = np.column_stack((id,id+1))

    equil_dist       = np.zeros(id.shape[0])
    equil_dist[0::3] = 1.453
    equil_dist[1::3] = 1.526
    equil_dist[2::3] = 1.300

    if n_chains > 1:
        between_chain = find_multichain_terms(id, chain_starts)
        ndx = np.where(between_chain)[0]
        id = id[ndx]
        equil_dist = equil_dist[ndx]
    create_array(dgrp, 'id', obj=id)

    grp = t.create_group(potential, 'Spring_bond')
    grp._v_attrs.arguments = np.array([b'Distance3D'])
    grp._v_attrs.integrator_level = 0
    grp._v_attrs.dim1 = 0
    grp._v_attrs.pbc  = 0
    grp._v_attrs.box_len  = 0

    spring_const = args.bond_stiffness*np.ones(id.shape[0])
    create_array(grp, 'id',           obj=np.arange(id.shape[0]))
    create_array(grp, 'equil_dist',   obj=equil_dist)
    create_array(grp, 'spring_const', obj=spring_const)

def write_angle_spring(args):

    agrp = t.create_group(potential, 'Angle')
    agrp._v_attrs.arguments = np.array([b'pos'])
    id = np.arange(n_atom-2)
    id = np.column_stack((id,id+2,id+1))

    equil_angles       = np.zeros(id.shape[0])
    equil_angles[0::3] = np.cos(109.5*deg)  # N->CA->C angle
    equil_angles[1::3] = np.cos(120.0*deg)  # CA->C->N angle
    equil_angles[2::3] = np.cos(120.0*deg)  # C->N->CA angle

    if n_chains > 1:
        between_chain = find_multichain_terms(id, chain_starts)
        ndx = np.where(between_chain)[0]
        id  = id[ndx]
        equil_angles = equil_angles[ndx]
    create_array(agrp, 'id', obj=id)

    grp = t.create_group(potential, 'Spring_angle')
    grp._v_attrs.arguments = np.array([b'Angle'])
    grp._v_attrs.integrator_level = 0
    grp._v_attrs.dim1 = 0
    grp._v_attrs.pbc  = 0
    grp._v_attrs.box_len  = 0

    spring_const = args.angle_stiffness*np.ones(id.shape[0])

    create_array(grp, 'id',           obj=np.arange(id.shape[0]))
    create_array(grp, 'equil_dist',   obj=equil_angles)
    create_array(grp, 'spring_const', obj=spring_const)

def write_angle_spring2(args):

    agrp = t.create_group(potential, 'VectorAngle')
    agrp._v_attrs.arguments = np.array([b'pos'])
    id = np.arange(n_atom-2)
    id = np.column_stack((id,id+2,id+1))
    if n_chains > 1:
        between_chain = find_multichain_terms(id, chain_starts)
        ndx = np.where(between_chain)[0]
        id  = id[ndx]
    create_array(agrp, 'id', obj=id)

    grp = t.create_group(potential, 'Spring_angle')
    grp._v_attrs.arguments = np.array([b'VectorAngle'])
    grp._v_attrs.integrator_level = 0
    grp._v_attrs.dim1 = 0
    grp._v_attrs.pbc  = 0
    grp._v_attrs.box_len  = 0

    spring_const = args.angle_stiffness*np.ones(id.shape[0])
    equil_angles = np.zeros(id.shape[0])
    equil_angles[0::3] = np.cos(109.5*deg)  # N->CA->C angle
    equil_angles[1::3] = np.cos(120.0*deg)  # CA->C->N angle
    equil_angles[2::3] = np.cos(120.0*deg)  # C->N->CA angle

    create_array(grp, 'id',           obj=np.arange(id.shape[0]))
    create_array(grp, 'equil_dist',   obj=equil_angles)
    create_array(grp, 'spring_const', obj=spring_const)

def write_omega_spring1(args, fasta_seq):
    # this is primarily used for omega bonds
    dgrp = t.create_group(potential, 'Dihedral_omega')
    dgrp._v_attrs.arguments = np.array([b'pos'])

    id = np.arange(1,n_atom-3,3)  # start at CA atom
    id = np.column_stack((id,id+1,id+2,id+3))
    ndx = np.arange(id.shape[0])
    if n_chains > 1:
        between_chain = find_multichain_terms(id, chain_starts)
        ndx = np.where(between_chain)[0]
        id  = id[ndx]
    create_array(dgrp, 'id', obj=id)

    grp = t.create_group(potential, 'Spring_omega')
    grp._v_attrs.arguments = np.array([b'Dihedral_omega'])
    grp._v_attrs.integrator_level = 0
    grp._v_attrs.dim1 = 0
    grp._v_attrs.pbc  = 1
    grp._v_attrs.box_len  = np.pi*2
    target_angle = np.where((fasta_seq[1:]=='CPR'), 0.*deg, 180.*deg)
    create_array(grp, 'id',           obj=np.arange(id.shape[0]))
    create_array(grp, 'equil_dist',   obj=target_angle[ndx])
    create_array(grp, 'spring_const', obj=args.omega_stiffness*np.ones(id.shape[0]))

def write_omega_spring2(parser, args, fasta_seq, pro_state_file):

    n_res = len(fasta_seq)

    fields = [ln.split() for ln in open(pro_state_file)]

    header_fields = 'residue init_state energy k width'.split()
    if [x.lower() for x in fields[0]] == header_fields:
        len_header_fields = len(header_fields)
    else:
        parser.error('First line of restraint table must be "%s"'%(" ".join(header_fields)))
    if not all(len(f)==len_header_fields for f in fields):
        parser.error('Invalid format for restraint file')

    fields = fields[1:]
    n_pro = len(fields)
    
    pro_res_id = np.zeros(n_pro, dtype=int)
    pro_state  = np.zeros(n_pro, dtype=int)
    energy     = np.zeros(n_pro)
    k          = np.zeros(n_pro)
    width      = np.zeros(n_pro)

    for i,f in enumerate(fields):
        pro_res_id[i] = int(f[0])
        msg = 'Double wall energy specified for residue %i (zero is first residue) but there are only %i residues in the FASTA'
        if not (0 <= pro_res_id[i] < n_res): raise ValueError(msg % (pro_res_id[i], n_res))
        pro_state[i]  = int(f[1])
        energy[i]     = float(f[2])  # the height from the trans state to the transition state
        k[i]          = float(f[3])  # k*energy is the  height from the cis state to the transition state
        width[i]      = float(f[4])  # compact_sigmoid cuts off at angle +/- width (degree)
        if width[i] <= 0.: raise ValueError('Cannot have negative width value')

    highlight_residues('residues that participate in any double wall potential in uppercase', fasta_seq, pro_res_id.ravel())
    # 0-based indexing sometimes trips up users, so give them a quick check
    S = 1/(width*deg)

    ndx1 = pro_res_id*1
    ndx0 = []
    for n in range(n_res):
        if not (n in ndx1):
            ndx0.append(n)
    ndx0 = np.array(ndx0)

    # this is primarily used for omega bonds
    id = np.arange(1,n_atom-3,3)  # start at CA atom
    id = np.column_stack((id,id+1,id+2,id+3))
    if n_chains > 1:
        between_chain = find_multichain_terms(id, chain_starts)
        ndx = np.where(between_chain)[0]
        id  = id[ndx]

    id0 = []
    id1 = []
    target_angle = []
    for i,ii in enumerate(id):
        rid = (ii[3]-1)//3
        if rid in ndx0:
            id0.append(id[i])
            if fasta_seq[rid]=='CPR':
                 target_angle.append(0.*deg)
            else:
                 target_angle.append(180.*deg)
        if rid in ndx1:
            id1.append(id[i])
    id0 = np.array(id0)
    id1 = np.array(id1)
    target_angle = np.array(target_angle)

    dgrp = t.create_group(potential, 'Dihedral_omega')
    dgrp._v_attrs.arguments = np.array([b'pos'])
    create_array(dgrp, 'id', obj=id0)

    grp = t.create_group(potential, 'Spring_omega')
    grp._v_attrs.arguments = np.array([b'Dihedral_omega'])
    grp._v_attrs.integrator_level = 0
    grp._v_attrs.dim1 = 0
    grp._v_attrs.pbc  = 1
    grp._v_attrs.box_len  = np.pi*2
    create_array(grp, 'id',           obj=np.arange(id0.shape[0]))
    create_array(grp, 'equil_dist',   obj=target_angle)
    create_array(grp, 'spring_const', obj=args.omega_stiffness*np.ones(id0.shape[0]))

    dgrp = t.create_group(potential, 'Dihedral_OmegaTransCis')
    dgrp._v_attrs.arguments = np.array([b'pos'])
    create_array(dgrp, 'id', obj=id1)

    one = np.ones(n_pro)
    ids = np.arange(n_pro)
    border = np.linspace(-np.pi, np.pi, 9)
    a = [ energy, energy*(1.-k),    energy, energy*0. ]
    b = [-energy,      energy*k, -energy*k, energy    ]
    node_name = ['trans1', 'cis0', 'cis1', 'trans0']
    for si in range(4):
        grp = t.create_group(potential, 'SigmoidEnergy_{}'.format(node_name[si]))
        grp._v_attrs.arguments = np.array([b'Dihedral_OmegaTransCis'])
        grp._v_attrs.integrator_level = 0
        grp._v_attrs.use_cutoff = 1
        grp._v_attrs.dim1 = 0
        create_array( grp, 'id',  obj = ids                )
        create_array( grp, 'a',   obj = a[si]              )
        create_array( grp, 'b',   obj = b[si]              )
        create_array( grp, 'c',   obj = one*border[si*2+1] )
        create_array( grp, 's',   obj = S                  )
        create_array( grp, 'min', obj = one*border[si*2]   )
        create_array( grp, 'max', obj = one*border[si*2+2] )

#---------------------------------------------------------------------------
#                               RAMA
#---------------------------------------------------------------------------

def mixture_potential(weights, potentials):
    ''' potentials must be normalized to the same value, preferably 1 '''
    potentials = np.array(potentials)
    assert len(weights) == len(potentials)
    weights = np.array(weights)
    weights = weights / weights.sum(axis=0)

    # ensure that we broadcast correctly against the potential
    weight_broadcast_shape = weights.shape + (1,)*(len(potentials.shape)-len(weights.shape))
    weights = weights.reshape(weight_broadcast_shape)

    potentials = potentials - np.log(weights)

    min_pot = potentials.min(axis=0)
    return min_pot - np.log(np.exp(min_pot-potentials).sum(axis=0))

def read_rama_maps_and_weights(seq, rama_group, mode='mixture', allow_CPR=True):
    assert mode in ['mixture', 'product']
    restype = rama_group._v_attrs.restype
    dirtype = rama_group._v_attrs.dir

    if type(restype[0]) == str:
        ridx_dict = dict([(x,i) for i,x in enumerate(restype)])
    else:
        ridx_dict = dict([(x.decode('utf-8'),i) for i,x in enumerate(restype)])
    if type(restype[0]) == str:
        didx = dict([(x,i) for i,x in enumerate(dirtype)])
    else:
        didx = dict([(x.decode('utf-8'),i) for i,x in enumerate(dirtype)])
    ridx = lambda resname, keep_cpr=True: (ridx_dict[resname] if resname!='CPR' or keep_cpr else ridx_dict['PRO'])

    dimer_pot    = rama_group.dimer_pot[:]
    dimer_weight = rama_group.dimer_weight[:]

    assert len(seq) >= 3   # avoid bugs

    # cis-proline is only CPR when it is the central residue, otherwise just use PRO

    V = lambda r,d,n: dimer_pot   [ridx(r,allow_CPR), didx[d], ridx(n,False)]
    W = lambda r,d,n: dimer_weight[ridx(r,allow_CPR), didx[d], ridx(n,False)]

    pots    = np.zeros((len(seq), dimer_pot.shape[-2], dimer_pot.shape[-1]), dtype='f4')
    weights = np.zeros((len(seq),),dtype='f4')

    pots   [0] = V(seq[0], 'right', seq[1])
    weights[0] = W(seq[0], 'right', seq[1])

    for i,l,c,r in zip(range(1,len(seq)-1), seq[:-2], seq[1:-1], seq[2:]):
        if   mode == 'product':
            pots[i]    = V(c,'left',l) + V(c,'right',r) - V(c,'right','ALL')
            weights[i] = 0.5*(W(c,'left',l) + W(c,'right',r))  # always just average weights
        elif mode == 'mixture':
            # it's a little sticky to figure out what the mixing proportions should be
            # there is basically a one-sided vs two-sided problem (what if we terminate a sheet?)
            # I am going with one interpretation that may not be right
            pots[i]    = mixture_potential([W(c,'left',l), W(c,'right',r)], [V(c,'left',l), V(c,'right',r)])
            weights[i] = 0.5*(W(c,'left',l) + W(c,'right',r))
        else:
            raise RuntimeError('impossible')

    pots   [-1] = V(seq[-1], 'left', seq[-2])
    weights[-1] = W(seq[-1], 'left', seq[-2])

    # Ensure normalization
    pots -= -np.log(np.exp(-1.0*pots).sum(axis=(-2,-1), keepdims=1))

    return pots, weights

def read_weighted_maps(seq, rama_library_h5, sheet_mixing=None, mode='mixture'):
    with tb.open_file(rama_library_h5) as tr:
        coil_pots, coil_weights = read_rama_maps_and_weights(seq, tr.root.coil, mode=mode)

        if sheet_mixing is None:
            return coil_pots
        else:
            sheet_pots, sheet_weights = read_rama_maps_and_weights(seq, tr.root.sheet, allow_CPR=False)
            return mixture_potential([coil_weights, sheet_weights*np.exp(-sheet_mixing)],
                                     [coil_pots,    sheet_pots])

def write_rama_map_pot(seq, rama_library_h5, sheet_mixing_energy=None, secstr_bias='', mode='mixture', param_deriv=False):
    grp = t.create_group(potential, 'rama_map_pot')
    grp._v_attrs.arguments = np.array([b'rama_coord'])
    grp._v_attrs.integrator_level = 0

    # create a sheet array for all the residures
    with tb.open_file(rama_library_h5) as tr:
        sheet_restype    = tr.root.sheet._v_attrs.restype
        sheet_rid_dict   = dict([(x,i) for i,x in enumerate(sheet_restype)])
        grp._v_attrs.restype = sheet_restype

    if sheet_mixing_energy is not None:
        # support finite differencing for potential derivative
        eps = 1e-2
        grp._v_attrs.sheet_eps = eps
        sheet_mixing_values = np.loadtxt(sheet_mixing_energy)
        assert sheet_mixing_values.size == len(sheet_restype)

    sheet      = []
    sheet_rids = []

    for s in seq:
        if s == "CPR":
            rid = sheet_rid_dict["PRO"]
        else:
            rid = sheet_rid_dict[s]
        sheet_rids.append(rid)
        sheet.append(sheet_mixing_values[rid])

    sheet      = np.array(sheet)
    sheet_rids = np.array(sheet_rids)
    rama_pot   = read_weighted_maps(seq, rama_library_h5, sheet, mode)

    # support finite differencing for potential derivative
    if param_deriv:
        eps = 5e-4
        grp._v_attrs.sheet_eps = eps
        for s in sheet_restype:

            if s == "PRO":
                if (s not in seq) and ("CPR" not in seq):
                    continue
            else:
                if s not in seq:
                    continue

            rid = sheet_rid_dict[s]

            more_sheet = sheet.copy()
            more_sheet[sheet_rids == rid] += eps
            less_sheet = sheet.copy()
            less_sheet[sheet_rids == rid] -= eps
            create_array(grp, 'more_sheet_rama_pot_'+s, read_weighted_maps(seq, rama_library_h5, more_sheet))
            create_array(grp, 'less_sheet_rama_pot_'+s, read_weighted_maps(seq, rama_library_h5, less_sheet))
 
    if secstr_bias:
        assert len(rama_pot.shape) == 3
        phi = np.linspace(-np.pi,np.pi,rama_pot.shape[1],endpoint=False)[:,None]
        psi = np.linspace(-np.pi,np.pi,rama_pot.shape[2],endpoint=False)[None,:]
        sigmoid_lessthan = lambda a,b: 1./(1.+np.exp(-(b-a)/(10.*deg)))

        helical_basin = sigmoid_lessthan(phi,0.*deg) * sigmoid_lessthan(-100.*deg,psi) * sigmoid_lessthan(psi,50.*deg)
        sheet_basin   = sigmoid_lessthan(phi,0.*deg) * (sigmoid_lessthan(psi,-100.*deg) + sigmoid_lessthan(50.*deg,psi))

        f = (ln.split() for ln in open(secstr_bias))
        assert f.next() == 'residue secstr energy'.split()
        for residue,secstr,energy in f:
            residue = int(residue)
            energy = float(energy)

            if secstr == 'helix':
                rama_pot[residue] += energy * helical_basin
            elif secstr == 'sheet':
                rama_pot[residue] += energy *   sheet_basin
            else:
                raise ValueError('secstr in secstr-bias file must be helix or sheet')

    # let's remove the average energy from each Rama map
    # so that the Rama potential emphasizes its variation

    rama_pot -= (rama_pot*np.exp(-rama_pot)).sum(axis=(-2,-1),keepdims=1)

    create_array(grp, 'residue_id',   obj=np.arange(len(seq)))
    create_array(grp, 'rama_map_id',  obj=np.arange(rama_pot.shape[0]))
    create_array(grp, 'rama_map_id_all',  obj=np.arange(rama_pot.shape[0]))
    create_array(grp, 'rama_pot',     obj=rama_pot)

def write_rama_map_pot2(parser, seq, rama_library_h5, pro_state_file, sheet_mixing_energy=None, mode='mixture' ):

    n_res = len(seq)
    seq_new = np.array(seq)

    fields = [ln.split() for ln in open(pro_state_file)]

    header_fields = 'residue init_state energy k width'.split()
    if [x.lower() for x in fields[0]] == header_fields:
        len_header_fields = len(header_fields)
    else:
        parser.error('First line of restraint table must be "%s"'%(" ".join(header_fields)))
    if not all(len(f)==len_header_fields for f in fields):
        parser.error('Invalid format for restraint file')

    fields = fields[1:]
    n_pro = len(fields)
    
    pro_res_id = np.zeros(n_pro, dtype=int)
    pro_state  = np.zeros(n_pro, dtype=int)
    for i,f in enumerate(fields):
        pro_res_id[i] = int(f[0])
        pro_state[i]  = int(f[1])

    state_list = ["PRO", "CPR"]
    for i,pr in enumerate(pro_res_id):
        is_pro = (seq[pr] in state_list)
        assert(is_pro)
        #seq_new[pr] = state_list[pro_state[i]]
        seq_new[pr] = state_list[0]

    # create a sheet array for all the residures
    with tb.open_file(rama_library_h5) as tr:
        sheet_restype  = tr.root.sheet._v_attrs.restype
        sheet_rid_dict = dict([(x,i) for i,x in enumerate(sheet_restype)])

    if sheet_mixing_energy is not None:
        # support finite differencing for potential derivative
        sheet_mixing_values = np.loadtxt(sheet_mixing_energy)
        assert sheet_mixing_values.size == len(sheet_restype)

    sheet = []
    for s in seq_new:
        if s == "CPR":
            rid = sheet_rid_dict["PRO"]
        else:
            rid = sheet_rid_dict[s]
        sheet.append(sheet_mixing_values[rid])
    sheet = np.array(sheet)

    rama_pot   = read_weighted_maps(seq_new, rama_library_h5, sheet, mode)
    rama_pot  -= (rama_pot*np.exp(-rama_pot)).sum(axis=(-2,-1),keepdims=1)

    seq_new2 = seq_new[:]
    for i,r in enumerate(pro_res_id):
        seq_new2[r] = state_list[1]

    ndx1 = pro_res_id*1
    ndx0 = []
    for n in range(n_res):
        if not (n in ndx1):
            ndx0.append(n)
    ndx0 = np.array(ndx0)

    grp = t.create_group( potential, 'rama_map_pot' )
    grp._v_attrs.arguments = np.array([b'rama_coord'])
    grp._v_attrs.integrator_level = 0
    grp._v_attrs.restype = sheet_restype
    create_array(grp, 'residue_id',      obj=np.arange(n_res)[ndx0])
    create_array(grp, 'rama_map_id',     obj=np.arange(rama_pot.shape[0])[ndx0])
    create_array(grp, 'rama_map_id_all', obj=np.arange(rama_pot.shape[0]))
    create_array(grp, 'rama_pot',        obj=rama_pot)

    rama_pot_trans  = read_weighted_maps(seq_new, rama_library_h5, sheet, mode)
    rama_pot_trans -= (rama_pot_trans*np.exp(-rama_pot_trans)).sum(axis=(-2,-1),keepdims=1)
    rama_pot_cis    = read_weighted_maps(seq_new2, rama_library_h5, sheet, mode)
    rama_pot_cis   -= (rama_pot_cis*np.exp(-rama_pot_cis)).sum(axis=(-2,-1),keepdims=1)

    grp = t.create_group(potential, 'SigmoidCoord_trans1')
    grp._v_attrs.arguments = np.array([b'Dihedral_OmegaTransCis'])
    grp._v_attrs.dim1 = 0
    create_array(grp, 'id',   obj=np.arange(ndx1.size))
    create_array(grp, 'c',    obj=np.ones(ndx1.size)*-0.5*np.pi)
    create_array(grp, 's',    obj=np.ones(ndx1.size)*3.0/np.pi)
    create_array(grp, 'sign', obj=np.ones(ndx1.size))

    grp = t.create_group(potential, 'SigmoidCoord_trans2')
    grp._v_attrs.arguments = np.array([b'Dihedral_OmegaTransCis'])
    grp._v_attrs.dim1 = 0
    create_array(grp, 'id',   obj=np.arange(ndx1.size))
    create_array(grp, 'c',    obj=np.ones(ndx1.size)*0.5*np.pi)
    create_array(grp, 's',    obj=np.ones(ndx1.size)*3.0/np.pi)
    create_array(grp, 'sign', obj=np.ones(ndx1.size)*-1.)

    grp = t.create_group(potential, 'Add_lambda_trans')
    grp._v_attrs.arguments = np.array([b'SigmoidCoord_trans1', b'SigmoidCoord_trans2'])
    grp._v_attrs.n_dim = 1
    grp._v_attrs.n_size = ndx1.size
    create_array(grp, 'dim1',    obj=np.array([0]))
    create_array(grp, 'dim2',    obj=np.array([0]))
    create_array(grp, 'id_pos1', obj=np.arange(ndx1.size))
    create_array(grp, 'id_pos2', obj=np.arange(ndx1.size))
    create_array(grp, 'id1_out', obj=np.arange(ndx1.size))
    create_array(grp, 'id2_out', obj=np.arange(ndx1.size))

    grp = t.create_group(potential, 'SigmoidCoord_cis1')
    grp._v_attrs.arguments = np.array([b'Dihedral_OmegaTransCis'])
    grp._v_attrs.dim1 = 0
    create_array(grp, 'id',   obj=np.arange(ndx1.size))
    create_array(grp, 'c',    obj=np.ones(ndx1.size)*-0.5*np.pi)
    create_array(grp, 's',    obj=np.ones(ndx1.size)*3.0/np.pi)
    create_array(grp, 'sign', obj=np.ones(ndx1.size)*-1.)

    grp = t.create_group(potential, 'SigmoidCoord_cis2')
    grp._v_attrs.arguments = np.array([b'Dihedral_OmegaTransCis'])
    grp._v_attrs.dim1 = 0
    create_array(grp, 'id',   obj=np.arange(ndx1.size))
    create_array(grp, 'c',    obj=np.ones(ndx1.size)*0.5*np.pi)
    create_array(grp, 's',    obj=np.ones(ndx1.size)*3.0/np.pi)
    create_array(grp, 'sign', obj=np.ones(ndx1.size))

    grp = t.create_group(potential, 'Multiply_lambda_cis')
    grp._v_attrs.arguments = np.array([b'SigmoidCoord_cis1', b'SigmoidCoord_cis2'])
    grp._v_attrs.n_dim = 1
    grp._v_attrs.n_size = ndx1.size
    create_array(grp, 'dim1',    obj=np.array([0]))
    create_array(grp, 'dim2',    obj=np.array([0]))
    create_array(grp, 'id_pos1', obj=np.arange(ndx1.size))
    create_array(grp, 'id_pos2', obj=np.arange(ndx1.size))
    create_array(grp, 'id1_out', obj=np.arange(ndx1.size))
    create_array(grp, 'id2_out', obj=np.arange(ndx1.size))

    grp = t.create_group(potential, 'RamaMap2_trans')
    grp._v_attrs.arguments = np.array([b'rama_coord', b'Add_lambda_trans'])
    grp._v_attrs.integrator_level = 0
    create_array(grp, 'residue_id',   obj=np.arange(n_res)[ndx1])
    create_array(grp, 'rama_map_id',  obj=np.arange(n_res)[ndx1])
    create_array(grp, 'rama_pot',     obj=rama_pot_trans)

    grp = t.create_group(potential, 'RamaMap2_cis')
    grp._v_attrs.arguments = np.array([b'rama_coord', b'Multiply_lambda_cis'])
    grp._v_attrs.integrator_level = 0
    create_array(grp, 'residue_id',   obj=np.arange(n_res)[ndx1])
    create_array(grp, 'rama_map_id',  obj=np.arange(n_res)[ndx1])
    create_array(grp, 'rama_pot',     obj=rama_pot_cis)

#---------------------------------------------------------------------------
#                             pair interactions
#---------------------------------------------------------------------------

def write_backbone_pair(fasta):
    n_res = len(fasta)
    grp = t.create_group(potential, 'backbone_pairs')
    grp._v_attrs.arguments = np.array([b'affine_alignment'])
    grp._v_attrs.integrator_level = 1

    ref_pos = np.zeros((n_res,4,3))
    ref_pos[:,0] = (-1.19280531, -0.83127186,  0.)        # N
    ref_pos[:,1] = ( 0.,          0.,          0.)        # CA
    ref_pos[:,2] = ( 1.25222632, -0.87268266,  0.)        # C
    ref_pos[:,3] = ( 0.,          0.94375626,  1.2068012) # CB
    ref_pos[fasta=='GLY',3] = np.nan

    ref_pos -= ref_pos[:,:3].mean(axis=1)[:,None]

    create_array(grp, 'id', obj=np.arange(n_res))
    create_array(grp, 'ref_pos', obj=ref_pos)
    create_array(grp, 'n_atom',  obj=np.isfinite(grp.ref_pos[:].sum(axis=-1)).sum(axis=-1))


def write_count_hbond(fasta, loose_hbond ):
    n_res = len(fasta)

    infer_group = t.get_node('/input/potential/infer_H_O')

    n_donor    = infer_group.donors   .id.shape[0]
    n_acceptor = infer_group.acceptors.id.shape[0]

    igrp = t.create_group(potential, 'protein_hbond')
    igrp._v_attrs.arguments = np.array([b'infer_H_O'])

    if use_intensive_memory:
        igrp._v_attrs.max_n_edge = n_res * 8

    # group1 is the HBond donors
    create_array(igrp, 'index1', np.arange(0,n_donor))
    create_array(igrp, 'type1',  np.zeros(n_donor, dtype='i'))
    create_array(igrp, 'id1',    infer_group.donors.residue[:])

    # group 2 is the HBond acceptors
    create_array(igrp, 'index2', np.arange(n_donor,n_donor+n_acceptor))
    create_array(igrp, 'type2',  np.zeros(n_acceptor, dtype='i'))
    create_array(igrp, 'id2',    infer_group.acceptors.residue[:])

    # parameters are inner_barrier, inner_scale, outer_barrier, outer_scale, wall_dp, inv_dp_width
    create_array(igrp, 'interaction_param', np.array([[
        [(0.5   if loose_hbond else 1.4  ), 1./0.10,
         (3.1   if loose_hbond else 2.5  ), 1./0.125,
         (0.182 if loose_hbond else 0.682), 1./0.05,
         0.,   0.]]]))

def write_short_hbond(fasta, hbond_energy):
    n_res = len(fasta)
    infer_group = t.get_node('/input/potential/infer_H_O')
    d_residues = infer_group.donors   .residue[:]
    a_residues = infer_group.acceptors.residue[:]

    grp = t.create_group(potential, 'hbond_energy')
    grp._v_attrs.arguments = np.array([b'protein_hbond', b'rama_coord'])
    grp._v_attrs.integrator_level = 1

    with tb.open_file(hbond_energy) as data:
        params = data.root.parameter[:]

    for hbe in params[:3]:
        if hbe > 0.:
            print ('\n**** WARNING ****  hydrogen bond formation energy set to repulsive value\n')
            break

    #sharpness = 1.0/(np.pi/12.)
    #params = [helix_hbond_energy, sheet_hbond_energy, turn_hbond_energy,
    #          0.0,          sharpness, np.pi*11./12., sharpness,  # for turn or sheet/helix
    #          -np.pi*2./3., sharpness, np.pi/3.,      sharpness ] # for sheet or helix

    create_array(grp, 'parameters',     params)
    create_array(grp, 'donor_resid',    d_residues)
    create_array(grp, 'acceptor_resid', a_residues)
    create_array(grp, 'rama_resid',     obj=np.arange(n_res))

def write_rotamer_placement(fasta, placement_library, dynamic_placement, dynamic_1body, fix_rotamer, excluded_residues):
    def compute_chi1_state(angles):
        chi1_state = np.ones(angles.shape, dtype='i')
        chi1_state[(   0.*deg<=angles)&(angles<120.*deg)] = 0
        chi1_state[(-120.*deg<=angles)&(angles<  0.*deg)] = 2
        return chi1_state

    with tb.open_file(placement_library) as data:
        restype_num = dict((aa.decode('ASCII'),i) for i,aa in enumerate(data.root.restype_order[:]))

        if dynamic_placement:
            placement_pos = data.root.rotamer_center[:].transpose((2,0,1,3)) # put layer index first
        else:
            placement_pos = data.root.rotamer_center_fixed[:]

        if dynamic_1body:
            placement_energy = -np.log(data.root.rotamer_prob[:].transpose((2,0,1)))[...,None]
        else:
            placement_energy = data.root.rotamer_prob_fixed[:][...,None]

        start_stop = data.root.rotamer_start_stop_bead[:]
        find_restype =                       data.root.restype_and_chi_and_state[:,0].astype('i')
        find_chi1 =                          data.root.restype_and_chi_and_state[:,1]
        find_chi1_state = compute_chi1_state(data.root.restype_and_chi_and_state[:,1])
        find_chi2 =                          data.root.restype_and_chi_and_state[:,2]
        find_state =                         data.root.restype_and_chi_and_state[:,3].astype('i')

    fix = dict()
    if fix_rotamer:
        fields = [x.split() for x in list(open(fix_rotamer))]

        header = 'residue restype chain resnum chi1 chi2'
        actual_header = [x.lower() for x in fields[0]]
        if actual_header != header.split():
            raise RuntimeError('First line of fix-rotamer table must be "%s" but is "%s" for file %s'
                    %(header," ".join(actual_header),fix_rotamer))

        for residue, restype, chain, resnum, chi1, chi2 in fields[1:]:
            if fasta[int(residue)] != (restype if restype != 'CPR' else 'PRO'):
                raise RuntimeError("fix-rotamer file does not match FASTA"
                    + ", residue %i should be %s but fix-rotamer file has %s"%(
                        int(residue), fasta[int(residue)], restype))
            chi1 = float(chi1)*deg  # convert to radians internally
            chi2 = float(chi2)*deg

            if restype == 'GLY' or restype == 'ALA':
                fix_state = 0
            else:
                if np.isnan(chi1): continue
                # determine states that have the right restype and compatible chi1
                chi1_state = compute_chi1_state(np.array([chi1]))[0]
                restype_admissible = find_restype == restype_num[fasta[int(residue)]]
                chi1_admissible = find_chi1_state == chi1_state
                admissible = restype_admissible&chi1_admissible
                admissible_chi2 = find_chi2[admissible]
                admissible_state = find_state[admissible]
                if len(admissible_state)==1:  # handle short residues (like VAL)
                    fix_state = admissible_state[0]
                else:
                    if np.isnan(chi2): continue
                    # now find the closest chi2 among those states and read off the state index
                    chi2_dist = (admissible_chi2-chi2)%(2*np.pi)
                    chi2_dist[chi2_dist>np.pi] -= 2*np.pi  # find closest periodic image
                    fix_state = admissible_state[np.argmin(chi2_dist)]

            fix[int(residue)] = fix_state

    rama_residue = []
    affine_residue = []
    layer_index = []
    beadtype_seq = []
    id_seq = []
    ref_chi1_state = []

    count_by_n_rot = dict()

    for rnum,aa in enumerate(fasta):
        if rnum in excluded_residues:
            continue

        restype = restype_num[aa]
        start,stop,n_bead = start_stop[restype]
        assert (stop-start)%n_bead == 0
        n_rot = (stop-start)//n_bead

        # if it should be fixed, then we must modify these answers to get a single rotamer
        if rnum in fix:
            if not (0 <= fix[rnum] < n_rot): raise ValueError('invalid fix rotamer state')
            start,stop = start+n_bead*fix[rnum], start+n_bead*(fix[rnum]+1)
            n_rot = 1

        if n_rot not in count_by_n_rot:
            count_by_n_rot[n_rot] = 0;

        base_id = (count_by_n_rot[n_rot]<<n_bit_rotamer) + n_rot
        count_by_n_rot[n_rot] += 1

        rama_residue  .extend([rnum]*(stop-start))
        affine_residue.extend([rnum]*(stop-start))
        layer_index   .extend(np.arange(start,stop))
        beadtype_seq  .extend(['%s_%i'%(aa,i) for i in range(n_bead)]*n_rot)
        id_seq        .extend(np.arange(stop-start)//n_bead + (base_id<<n_bit_rotamer))

    sc_node_name = 'placement%s_point_vector_only' % ('' if dynamic_placement else '_fixed')
    grp = t.create_group(potential, sc_node_name)
    grp._v_attrs.arguments = np.array([b'affine_alignment'] + ([b'rama_coord'] if dynamic_placement else []))
    create_array(grp, 'rama_residue',    rama_residue)
    create_array(grp, 'affine_residue',  affine_residue)
    create_array(grp, 'layer_index',     layer_index)
    create_array(grp, 'placement_data',  placement_pos[...,:6])
    create_array(grp, 'beadtype_seq',    beadtype_seq)
    create_array(grp, 'id_seq',          np.array(id_seq))
    create_array(grp, 'fix_rotamer',     np.array(sorted(fix.items())))
    # create_array(grp, 'ref_chi1_state',  np.array(ref_chi1_state))
    # create_array(grp, 'find_chi1',       find_chi1)

    pl_node_name = 'placement%s_scalar' % ('' if dynamic_1body else '_fixed')
    grp = t.create_group(potential, pl_node_name)
    grp._v_attrs.arguments = np.array([b'affine_alignment']+([b'rama_coord'] if dynamic_1body else []))
    create_array(grp, 'rama_residue',    rama_residue)
    create_array(grp, 'affine_residue',  affine_residue)
    create_array(grp, 'layer_index',     layer_index)
    create_array(grp, 'placement_data',  placement_energy)

    return sc_node_name, pl_node_name


def write_rotamer_backbone(fasta, coverage_library, sc_node_name):

    n_res          = len(fasta)

    pre_PRO        = True
    backbone_group = backbone_group4
    group_id       = set(val for val in backbone_group.values())
    num_group      = len(group_id)

    # for donor and acceptor
    infer_group = t.get_node('/input/potential/infer_H_O')
    n_donor    = infer_group.donors   .id.shape[0]
    n_acceptor = infer_group.acceptors.id.shape[0]
    n_hbond    = n_donor + n_acceptor
    d_residues = infer_group.donors   .residue[:]
    a_residues = infer_group.acceptors.residue[:]

    # for sidechains
    sc_node = t.get_node(t.root.input.potential, sc_node_name)
    rseq      = sc_node.beadtype_seq[:]
    sc_resnum = sc_node.affine_residue[:]
    sc_index  = np.arange(len(rseq))

    with tb.open_file(coverage_library) as data:
         bead_num = dict((k,i) for i,k in enumerate(data.root.bead_order[:]))
         sc_type  = np.array([bead_num[s] for s in rseq])
         coverage_interaction   = data.root.coverage_interaction[:]
         hydrophobe_placement   = data.root.hydrophobe_placement[:]
         hydrophobe_interaction = data.root.hydrophobe_interaction[:]

    # sc-hbond interaction
    hb_index = np.arange(n_hbond)
    cgrp = t.create_group(potential, 'hbond_coverage')
    cgrp._v_attrs.arguments = np.array([b'protein_hbond', bstring(sc_node_name)])
    if use_intensive_memory:
        cgrp._v_attrs.max_n_edge = n_res * 80

    create_array(cgrp, 'interaction_param', coverage_interaction)
    # create the hbond type 
    hb_type = []
    for resid in d_residues:
        res     = fasta[resid]
        bbt     = backbone_group[res]
        if resid < n_res-1 and pre_PRO:
            if fasta[resid+1] in ['PRO', 'CPR']:
                bbt = num_group
        hb_type.append(bbt*2+0)
    for resid in a_residues:
        res     = fasta[resid]
        bbt     = backbone_group[res]
        if resid < n_res-1 and pre_PRO:
            if fasta[resid+1] in ['PRO', 'CPR']:
                bbt = num_group
        hb_type.append(bbt*2+1)
    hb_type = np.array(hb_type)
    # group1 is the HBond partners
    create_array(cgrp, 'index1', hb_index)
    create_array(cgrp, 'type1',  hb_type)  # donor is 0, acceptor is 1
    create_array(cgrp, 'id1',    np.concatenate([d_residues, a_residues]))
    # group2 is the sc
    create_array(cgrp, 'index2', sc_index)
    create_array(cgrp, 'type2',  sc_type)
    create_array(cgrp, 'id2',    sc_resnum)

    # create the oriented backbone atoms
    bb_index  = np.arange(3*n_res)
    bb_resnum = bb_index/3
    grp = t.create_group(potential, 'placement_fixed_point_vector_scalar')
    grp._v_attrs.arguments = np.array([b'affine_alignment'])
    create_array(grp, 'affine_residue',  bb_resnum)
    create_array(grp, 'layer_index',     bb_index%3)
    create_array(grp, 'placement_data',  hydrophobe_placement)

    # sc-backbone interaction
    cgrp = t.create_group(potential, 'hbond_coverage_hydrophobe')
    cgrp._v_attrs.arguments = np.array([b'placement_fixed_point_vector_scalar', bstring(sc_node_name)])
    if use_intensive_memory:
        cgrp._v_attrs.max_n_edge = n_res * 120
    create_array(cgrp, 'interaction_param', hydrophobe_interaction)
    # for backbone type
    bb_type = []
    for resid, res in enumerate(fasta):
        for i in [0,1,2]:
            bbt = backbone_group[res]
            if resid < n_res-1 and pre_PRO:
                if fasta[resid+1] in ['PRO', 'CPR']:
                    bbt = num_group
            bb_type.append(bbt*3+i)
    bb_type = np.array(bb_type)
    # group1 is the backbone partners
    create_array(cgrp, 'index1', bb_index)
    create_array(cgrp, 'type1',  bb_type)
    create_array(cgrp, 'id1',    bb_resnum)
    # group 2 is the sidechains
    create_array(cgrp, 'index2', sc_index)
    create_array(cgrp, 'type2',  sc_type)
    create_array(cgrp, 'id2',    sc_resnum)

def write_rotamer(fasta, interaction_library, damping, sc_node_name, pl_node_name, suffix=''):

    n_res = len(fasta)

    g = t.create_group(t.root.input.potential, 'rotamer%s' % suffix)
    args = [bstring(sc_node_name), bstring(pl_node_name)]
    def arg_maybe(nm):
        if nm in t.root.input.potential: args.append(bstring(nm))
    arg_maybe('hbond_coverage')
    arg_maybe('hbond_coverage_hydrophobe')

    g._v_attrs.arguments = np.array(args)
    g._v_attrs.integrator_level = 1
    g._v_attrs.max_iter = 1000
    g._v_attrs.tol      = 1e-3
    g._v_attrs.damping  = damping
    g._v_attrs.iteration_chunk_size = 2

    pg = t.create_group(g, "pair_interaction")
    if use_intensive_memory:
        pg._v_attrs.max_n_edge = n_res * 200

    with tb.open_file(interaction_library) as data:
         create_array(pg, 'interaction_param', data.root.pair_interaction[:])
         bead_num = dict((k,i) for i,k in enumerate(data.root.bead_order[:]))
         # pg._v_attrs.energy_cap = data.root._v_attrs.energy_cap_1body
         # pg._v_attrs.energy_cap_width = data.root._v_attrs.energy_cap_width_1body

    sc_node = t.get_node(t.root.input.potential, sc_node_name)
    rseq = sc_node.beadtype_seq[:]
    create_array(pg, 'index', np.arange(len(rseq)))
    create_array(pg, 'type',  np.array([bead_num[s] for s in rseq]))
    create_array(pg, 'id',    sc_node.id_seq[:])


#---------------------------------------------------------------------------
#                           multibody interactions
#---------------------------------------------------------------------------

def write_environment(fasta, environment_library, sc_node_name, potential_type=0, weights_number = 20, vector_CA_CO = False, use_bb_bl_for_hb=False, excluded_residues=[], rot_exclude_residues=[]):

    n_res = len(fasta)

    with tb.open_file(environment_library) as lib:
        coverage_param = lib.root.coverage_param[:]
        # params are r0,r_sharpness, dot0, dot_sharpness
        weights           = lib.root.weights[:]
        restype_order = dict([(x.decode('ASCII'),i) for i,x in enumerate(lib.root.restype_order[:])])

    # Place CB
    pgrp = t.create_group(potential, 'placement_fixed_point_vector_only_CB')
    pgrp._v_attrs.arguments = np.array([b'affine_alignment'])
    ref_pos = np.zeros((4,3))
    ref_pos[0] = (-1.19280531, -0.83127186,  0.)        # N
    ref_pos[1] = ( 0.,          0.,          0.)        # CA
    ref_pos[2] = ( 1.25222632, -0.87268266,  0.)        # C
    ref_pos[3] = ( 0.,          0.94375626,  1.2068012) # CB
    ref_pos   -= ref_pos[:3].mean(axis=0,keepdims=1)

    placement_data = np.zeros((1,6))
    placement_data[0,0:3] = ref_pos[3]

    if vector_CA_CO:
        placement_data[0,3:6] = (ref_pos[3]-ref_pos[2])/vmag(ref_pos[3]-ref_pos[2])
    else:
        placement_data[0,3:6] = (ref_pos[3]-ref_pos[1])/vmag(ref_pos[3]-ref_pos[1])   

    create_array(pgrp, 'affine_residue',  np.arange(len(fasta)))
    create_array(pgrp, 'layer_index',     np.zeros(len(fasta),dtype='i'))
    create_array(pgrp, 'placement_data',  placement_data)

    # Bring position and probability together for the side chains
    sc_node = t.get_node(t.root.input.potential, sc_node_name)
    n_sc = sc_node.affine_residue.shape[0]

    # Compute SC coverage of the CB
    cgrp = t.create_group(potential, 'environment_coverage_sc')
    cgrp._v_attrs.arguments = np.array([b'placement_fixed_point_vector_only_CB',b'weighted_pos'])
    cgrp._v_attrs.num_aa_types = 20
    if use_intensive_memory:
        cgrp._v_attrs.max_n_edge = n_res * 160

    res_seq = np.array([resid for resid in range(len(fasta)) if resid not in excluded_residues]) 

    # group1 is the source CB
    create_array(cgrp, 'index1', res_seq)
    create_array(cgrp, 'type1',  np.array([restype_order[s] for resid, s in enumerate(fasta) if resid not in excluded_residues]))  # one type per CB type
    create_array(cgrp, 'id1',    res_seq)
    # group 2 is the weighted points to interact with
    create_array(cgrp, 'index2', np.arange(n_sc))
    create_array(cgrp, 'type2',  0*np.arange(n_sc))   # for now coverage is very simple, so no types on SC
    create_array(cgrp, 'id2',    sc_node.affine_residue[:])
    create_array(cgrp, 'interaction_param', coverage_param)
    create_array(cgrp, 'aa_types', np.array([restype_order[s] for resid, s in enumerate(fasta) if resid not in excluded_residues]))

    # Couple an energy to the coverage coordinates
    assert (potential_type == 0 or potential_type == 1)

    coupling_method = "nonlinear" if potential_type == 0 else "sigmoid"
    egrp = t.create_group(potential, '{}_coupling_environment'.format(coupling_method))
    egrp._v_attrs.arguments = np.array([b'environment_coverage_sc'])
    egrp._v_attrs.integrator_level = 1
    egrp._v_attrs.number_independent_weights = weights_number # 1, or 20, or 400

    create_array(egrp, 'coupling_types', np.array([restype_order[s] for resid, s in enumerate(fasta) if resid not in excluded_residues]))
    create_array(egrp, 'weights', np.array(weights))

    # spline-like coupling function
    if potential_type == 0:
        with tb.open_file(environment_library) as lib:
            energies          = lib.root.energies[:]
            energies_x_offset = lib.root.energies._v_attrs.offset
            energies_x_inv_dx = lib.root.energies._v_attrs.inv_dx
        create_array(egrp, 'coeff', energies)
        egrp.coeff._v_attrs.spline_offset = energies_x_offset
        egrp.coeff._v_attrs.spline_inv_dx = energies_x_inv_dx
    # sigmoid-like coupling function
    else:
        with tb.open_file(environment_library) as lib:
            scales    = lib.root.scale[:]
            centers   = lib.root.center[:]
            sharpness = lib.root.sharpness[:]
        create_array(egrp, 'scale', np.array(scales))
        create_array(egrp, 'center', np.array(centers))
        create_array(egrp, 'sharpness', np.array(sharpness))

def write_bb_environment(fasta, environment_library, sc_node_name, bb_env_fn, use_hbb=False):

    n_res = len(fasta)

    with tb.open_file(environment_library) as lib:
        coverage_param = lib.root.coverage_param[:]
        # params are r0,r_sharpness, dot0, dot_sharpness
        weights           = lib.root.weights[:]
        restype_order = dict([(x.decode('ASCII'),i) for i,x in enumerate(lib.root.restype_order[:])])

    # Bring position and probability together for the side chains
    sc_node = t.get_node(t.root.input.potential, sc_node_name)
    n_sc = sc_node.affine_residue.shape[0]

    #
    infer_group = t.get_node('/input/potential/infer_H_O')
    d_res_id    = infer_group.donors.residue[:]
    a_res_id    = infer_group.acceptors.residue[:]
    res_id      = np.concatenate([d_res_id, a_res_id])
    n_donor     = infer_group.donors   .id.shape[0]
    n_acceptor  = infer_group.acceptors.id.shape[0]

    # Compute SC coverage of the HN or OC
    cgrp = t.create_group(potential, 'environment_coverage_hb')
    cgrp._v_attrs.arguments = np.array([b'infer_H_O',b'weighted_pos'])
    cgrp._v_attrs.num_aa_types = 20
    if use_intensive_memory:
        cgrp._v_attrs.max_n_edge = n_res * 160

    coverage_param = np.array([[[6.0, 1.0, -0.1, 1.0]]])
    create_array(cgrp, 'interaction_param', coverage_param)
    create_array(cgrp, 'aa_types', [restype_order[fasta[i]] for i in res_id])
    # group1 is the H or O
    create_array(cgrp, 'index1', np.arange(n_donor+n_acceptor))
    create_array(cgrp, 'type1',  np.arange(n_donor+n_acceptor)*0)
    create_array(cgrp, 'id1',    res_id)
    # group 2 is the weighted points to interact with
    create_array(cgrp, 'index2', np.arange(n_sc))
    create_array(cgrp, 'type2',  0*np.arange(n_sc))   # for now coverage is very simple, so no types on SC
    create_array(cgrp, 'id2',    sc_node.affine_residue[:])

    # Compute acceptor coverage of the HN
    cgrp = t.create_group(potential, 'hb_environment_coverage_hn')
    cgrp._v_attrs.arguments = np.array([b'infer_H_O',b'infer_H_O'])
    if use_intensive_memory:
        cgrp._v_attrs.max_n_edge = n_res * 6

    coverage_param = np.array([[[3.5, 2.0, 0.0, 2.0]]])
    create_array(cgrp, 'interaction_param', coverage_param)
    # group1 is the H
    create_array(cgrp, 'index1', np.arange(n_donor))
    create_array(cgrp, 'type1',  np.arange(n_donor)*0)
    create_array(cgrp, 'id1',    np.array(d_res_id))
    # group 2 is the O
    create_array(cgrp, 'index2', np.arange(n_donor, n_donor+n_acceptor))
    create_array(cgrp, 'type2',  np.arange(n_acceptor)*0)
    create_array(cgrp, 'id2',    np.array(a_res_id))

    # Compute donor coverage of the OC
    cgrp = t.create_group(potential, 'hb_environment_coverage_oc')
    cgrp._v_attrs.arguments = np.array([b'infer_H_O',b'infer_H_O'])
    if use_intensive_memory:
        cgrp._v_attrs.max_n_edge = n_res * 6

    coverage_param = np.array([[[3.5, 2.0, 0.0, 2.0]]])
    create_array(cgrp, 'interaction_param', coverage_param)
    # group1 is the O
    create_array(cgrp, 'index1', np.arange(n_donor, n_donor+n_acceptor))
    create_array(cgrp, 'type1',  np.arange(n_acceptor)*0)
    create_array(cgrp, 'id1',    np.array(a_res_id))
    # group 2 is the H
    create_array(cgrp, 'index2', np.arange(n_donor))
    create_array(cgrp, 'type2',  np.arange(n_donor)*0)
    create_array(cgrp, 'id2',    np.array(d_res_id))

    # Compute donor coverage of the BB
    if use_hbb:
        cgrp = t.create_group(potential, 'hbbb_coverage')
        cgrp._v_attrs.arguments = np.array([b'infer_H_O',b'pos'])
        coverage_param = np.array([[[6.0, 1.0, -0.1, 1.0]]])
        if use_intensive_memory:
            cgrp._v_attrs.max_n_edge = n_res * 200

        create_array(cgrp, 'interaction_param', coverage_param)
        # group1 is the O
        create_array(cgrp, 'index1', np.arange(n_donor))
        create_array(cgrp, 'type1',  np.arange(n_donor)*0)
        create_array(cgrp, 'id1',    np.array(d_res_id))
        # group 2 is the H
        create_array(cgrp, 'index2', np.arange(n_res*3))
        create_array(cgrp, 'type2',  np.arange(n_res*3)*0)
        create_array(cgrp, 'id2',    np.arange(n_res*3)/3)

    # cat them 
    pgrp = t.create_group(potential, 'cat_pos_bb_coverage')
    pgrp._v_attrs.arguments = np.array([b'hb_environment_coverage_hn', b'hb_environment_coverage_oc'])
    pgrp._v_attrs.n_dim = 1
    create_array(pgrp, 'index_pos1', np.arange(n_donor))
    create_array(pgrp, 'index_pos2', np.arange(n_acceptor))

    #
    egrp = t.create_group(potential, 'bb_sigmoid_coupling_environment')
    egrp._v_attrs.arguments = np.array([b'environment_coverage_hb', b'cat_pos_bb_coverage'])
    egrp._v_attrs.integrator_level = 1
    egrp._v_attrs.num_aa_types = 20
    bb_env_param = np.loadtxt(bb_env_fn)
    egrp._v_attrs.scale = bb_env_param[0]
    egrp._v_attrs.center = bb_env_param[1]
    egrp._v_attrs.sharpness = bb_env_param[2]
    egrp._v_attrs.hbond_weight = bb_env_param[3]
    create_array(egrp, 'weights', np.array(weights))

#---------------------------------------------------------------------------
#                               membrane
#---------------------------------------------------------------------------

def write_surface_coord(fasta_seq, method, thickness, included_list=[], zbuffer=0.0 ):

    n_res = len(fasta_seq)

    # cat the bb and cb atoms together
    pgrp = t.create_group(potential, 'cat_pos_bb_cb')
    pgrp._v_attrs.arguments = np.array([b'pos', b'placement_fixed_point_only_CB'])
    pgrp._v_attrs.n_dim = 3
    index_bb = np.arange(n_res*3)
    index_cb = np.arange(n_res)
    create_array(pgrp, 'index_pos1', index_bb)
    create_array(pgrp, 'index_pos2', index_cb)

    # the auxiliary points are used to fill the surface
    auxiliary_point = [ [ 0.000,  1.000],
                        [ 0.000,  0.500],
                        [ 0.866,  0.500],
                        [ 0.433,  0.250],
                        [ 0.866, -0.500],
                        [ 0.433, -0.250],
                        [ 0.000, -1.000],
                        [ 0.000, -0.500],
                        [-0.866, -0.500],
                        [-0.433, -0.250],
                        [-0.866,  0.500],
                        [-0.433,  0.250],
                        [ 0.433,  0.750],
                        [ 0.866,  0.000],
                        [ 0.433, -0.750],
                        [-0.433, -0.750],
                        [-0.866,  0.000],
                        [-0.433,  0.750] ]

    # surface node
    grp = t.create_group(t.root.input.potential, 'surface')
    grp._v_attrs.arguments = np.array([b'cat_pos_bb_cb' ])
    grp._v_attrs.n_rotation = 20
    grp._v_attrs.half_thickness = thickness*0.5
    grp._v_attrs.zbuffer = zbuffer
    grp._v_attrs.exclusion = 0

    grp._v_attrs.len_gridx = 3.0
    grp._v_attrs.len_gridy = 3.0
    grp._v_attrs.n_gridz = 6
    grp._v_attrs.n_smallgridx = 1
    grp._v_attrs.n_smallgridy = 1
    grp._v_attrs.n_smallgridz = 2
    grp._v_attrs.atom_radius = 1.6

    if method not in [0,1]:
        print('Error')
    grp._v_attrs.method_type = method

    if len(included_list) > 0:
        grp._v_attrs.exclusion = 1
        create_array(grp, 'included_list', np.array(included_list))

    create_array( grp, 'index_cb',        index_cb+n_res*3 )
    create_array( grp, 'auxiliary_point', np.array(auxiliary_point) )

def write_membrane_potential(
        fasta_seq, membrane_potential_fpath, membrane_thickness, environment_potential, membrane_exposed_criterion, membrane_exclude_residues, hbond_exclude_residues):

    grp = t.create_group(t.root.input.potential, 'membrane_potential')
    grp._v_attrs.arguments = np.array([b'placement_fixed_point_only_CB', b'environment_coverage_sc', b'protein_hbond'])
    grp._v_attrs.integrator_level = 0

    with tb.open_file(environment_potential) as lib:
        weights         = lib.root.weights[:]

    with tb.open_file(membrane_potential_fpath) as lib:
        resnames        = lib.root.names[:]
        cb_energy       = lib.root.cb_energy[:]
        cb_z_min        = lib.root.cb_energy._v_attrs.z_min
        cb_z_max        = lib.root.cb_energy._v_attrs.z_max
        thickness       = lib.root.cb_energy._v_attrs.thickness
        uhb_energy      = lib.root.uhb_energy[:]
        uhb_z_min       = lib.root.uhb_energy._v_attrs.z_min
        uhb_z_max       = lib.root.uhb_energy._v_attrs.z_max
        cov_midpoint    = lib.root.cov_midpoint[:]
        cov_sharpness   = lib.root.cov_sharpness[:]

    if membrane_exposed_criterion:
        with tb.open_file(membrane_exposed_criterion) as lib:
            cov_midpoint    = lib.root.cov_midpoint[:]
            cov_sharpness   = lib.root.cov_sharpness[:]
    
    #<----- ----- ----- ----- donor/acceptor res ids ----- ----- ----- ----->#
    # Note: hbond_excluded_residues is the same as in the function write_infer_H_O.
    n_res                = len(fasta_seq)
    donor_residue_ids    = np.array([i for i in range(n_res) if i>0       and i not in hbond_exclude_residues and fasta_seq[i]!='PRO'])
    acceptor_residue_ids = np.array([i for i in range(n_res) if i<n_res-1 and i not in hbond_exclude_residues])

    #<----- ----- ----- ----- make energy splines ----- ----- ----- ----->#
    import scipy.interpolate
    def extrapolated_spline(x0, y0):
        spline = scipy.interpolate.InterpolatedUnivariateSpline(x0,y0)
        def f(x, spline=spline):
            return np.select(
                    [(x<x0[0]),              (x>x0[-1]),              np.ones_like(x,dtype='bool')],
                    [np.zeros_like(x)+y0[0], np.zeros_like(x)+y0[-1], spline(x)])
        return f

    cb_z_lib           = np.linspace(cb_z_min, cb_z_max, cb_energy.shape[-1])
    cb_energy_splines  = [extrapolated_spline(cb_z_lib, ene) for ene in cb_energy]

    uhb_z_lib          = np.linspace(uhb_z_min, uhb_z_max, uhb_energy.shape[-1])
    uhb_energy_splines = [extrapolated_spline(uhb_z_lib, ene) for ene in uhb_energy]

    #<----- ----- ----- ----- make energy splines ----- ----- ----- ----->#
    # This step is necessary in case the supplied membrane thickness is not eaual to the thickness in the membrane potential file.
    default_half_thickness = thickness/2.
    half_thickness         = membrane_thickness/2.
    z_                     = np.linspace(-half_thickness - 15., half_thickness + 15., int((membrane_thickness+30.)/0.25)+1)

    # ensure that the potential is continuous at 0
    # spline(z-(half_thickness-default_half_thickness)) may not equal to spline(z+(half_thickness-default_half_thickness))
    membrane_cb_energies = np.zeros((len(cb_energy_splines), len(z_)))
    for ispl, spline in enumerate(cb_energy_splines):
        if half_thickness < default_half_thickness:
            delta_t = default_half_thickness - half_thickness
            delta_s = spline(delta_t) - spline(-delta_t)
            membrane_cb_energies[ispl] = np.select([(z_ < 0), (z_ >= 0.)],
                                                   [spline(z_-delta_t) + 0.5*delta_s, spline(z_+delta_t) - 0.5*delta_s])
        elif half_thickness > default_half_thickness:
            delta_t = half_thickness - default_half_thickness
            membrane_cb_energies[ispl] = np.select([
                (z_ <  -delta_t),
                (z_ >= -delta_t) & (z_ <= delta_t),
                (z_ >   delta_t)],
                [spline(z_+delta_t), spline(0), spline(z_-delta_t)])
        else:
            membrane_cb_energies[ispl] = spline(z_)

    membrane_uhb_energies = np.zeros((len(uhb_energy_splines), len(z_)))
    for ispl, spline in enumerate(uhb_energy_splines):
        if half_thickness < default_half_thickness:
            delta_t = default_half_thickness - half_thickness
            delta_s = spline(delta_t) - spline(-delta_t)
            membrane_uhb_energies[ispl] = np.select([(z_ < 0), (z_ >= 0.)],
                                                    [spline(z_-delta_t) + 0.5*delta_s, spline(z_+delta_t) - 0.5*delta_s])
        elif half_thickness > default_half_thickness:
            delta_t = half_thickness - default_half_thickness
            membrane_uhb_energies[ispl] = np.select([
                (z_ <  -delta_t),
                (z_ >= -delta_t) & (z_ <= delta_t),
                (z_ >   delta_t)],
                [spline(z_+delta_t), spline(0), spline(z_-delta_t)])
        else:
            membrane_uhb_energies[ispl] = spline(z_)

    #<----- ----- ----- ----- cb energy indices ----- ----- ----- ----->#
    # Note: there's a residue type, NON, in resnames for those excluded from membrane potential.
    # And there's a potential profile in cb_energy for NON, which is all zeros. 
    if set(membrane_exclude_residues).difference(range(len(fasta_seq))) != set():
        raise ValueError('Residue number', set(membrane_exclude_residues).difference(range(len(fasta_seq))), 'not valid')
    highlight_residues('membrane_exclude_residues', fasta_seq, membrane_exclude_residues)

    sequence = list(fasta_seq)
    for num in membrane_exclude_residues:
        sequence[num] = 'NON' 
    sequence = np.array(sequence)

    resname_to_num  = dict([(aa,i) for i,aa in enumerate(resnames)])
    residue_id      = np.array([i for i,aa in enumerate(sequence)])
    cb_energy_index = np.array([resname_to_num[aa] for aa in sequence])

    #<----- ----- ----- ----- write to grp ----- ----- ----- ----->#
    create_array(grp,             'cb_index', residue_id)
    create_array(grp,            'env_index', residue_id)
    create_array(grp,         'residue_type', cb_energy_index)
    create_array(grp,         'cov_midpoint', cov_midpoint)
    create_array(grp,        'cov_sharpness', cov_sharpness)
    create_array(grp,            'cb_energy', membrane_cb_energies)
    create_array(grp,              'weights', weights[:20])
    create_array(grp,           'uhb_energy', membrane_uhb_energies)
    create_array(grp,    'donor_residue_ids', donor_residue_ids)
    create_array(grp, 'acceptor_residue_ids', acceptor_residue_ids)
    grp. cb_energy._v_attrs.z_min = z_[ 0]
    grp. cb_energy._v_attrs.z_max = z_[-1]
    grp.uhb_energy._v_attrs.z_min = z_[ 0]
    grp.uhb_energy._v_attrs.z_max = z_[-1]

def write_membrane_potential3(fasta_seq, membrane_potential_fpath, membrane_thickness, membrane_exclude_residues, hbond_exclude_residues, use_curvature=False, curvature_radius=200, curvature_sign=1):

    half_thickness = membrane_thickness * 0.5

    with tb.open_file(membrane_potential_fpath) as lib:
        resnames       = lib.root.names[:]
        resname_to_num = dict([(aa.decode('ASCII'),i) for i,aa in enumerate(resnames)])
        cb_energy      = lib.root.cb_energy[:]
        hb_energy      = lib.root.hb_energy[:]
        cb_z_min       = lib.root._v_attrs.cb_z_min
        cb_z_max       = lib.root._v_attrs.cb_z_max
        hb_z_min       = lib.root._v_attrs.hb_z_min
        hb_z_max       = lib.root._v_attrs.hb_z_max

        cb_n_node      = cb_energy.shape[2]
        hb_n_node      = hb_energy.shape[2]

        n_bl           = cb_energy.shape[1]
        burial_nodes   = lib.root.burial_nodes[:]
        assert (burial_nodes.size+1 == n_bl)

        cb_dz = (cb_z_max-cb_z_min)/(cb_n_node-2)
        hb_dz = (hb_z_max-hb_z_min)/(hb_n_node-2)

        cb_start = cb_z_min - half_thickness
        hb_start = hb_z_min - half_thickness

    curvature_center = Const3DCoord1(t, 'curvature_center', 'pos', 0, 1, 2, [[0., 0., curvature_radius*-1*curvature_sign ]])

    # memb potential for CB atoms
    grp = t.create_group(t.root.input.potential, 'cb_membrane_potential')
    grp._v_attrs.arguments            = np.array([b'placement_fixed_point_only_CB', b'environment_coverage_sc', bstring(curvature_center)])
    grp._v_attrs.integrator_level     = 0
    grp._v_attrs.z_start              = cb_start
    grp._v_attrs.z_scale              = 1./cb_dz
    grp._v_attrs.bl_range_sharpness   = 1.0
    grp._v_attrs.left_right_node      = 0.0
    grp._v_attrs.left_right_sharpness = 1./cb_dz
    grp._v_attrs.half_thickness       = half_thickness

    grp._v_attrs.use_curvature        = int(use_curvature)
    grp._v_attrs.curvature_radius     = curvature_radius
    grp._v_attrs.curvature_sign       = curvature_sign

    print(grp._v_attrs.use_curvature, grp._v_attrs.curvature_radius, grp._v_attrs.curvature_sign)

    n_res = len(fasta_seq)
    used_res = []

    if len(membrane_exclude_residues) > 0:
        for r in range(n_res):
            if not (r in membrane_exclude_residues):
                used_res.append(r)
        used_res = np.array(used_res)
    else:
        used_res = np.arange(n_res)

    res_type = []
    for r in used_res:
        res_type.append(resname_to_num[fasta_seq[r]])
    res_type = np.array(res_type)

    env_group  = t.get_node('/input/potential/sigmoid_coupling_environment')
    weights    = env_group.weights[:20]

    create_array(grp,  'cb_index', used_res)
    create_array(grp,  'res_type', res_type)
    create_array(grp,   'weights', weights)
    create_array(grp,  'bl_nodes', burial_nodes)
    create_array(grp,     'coeff', cb_energy)

    # memb potential for HB atoms
    grp = t.create_group(t.root.input.potential, 'hb_membrane_potential')
    grp._v_attrs.arguments            = np.array([b'protein_hbond', bstring(curvature_center)])
    grp._v_attrs.integrator_level     = 0
    grp._v_attrs.z_start              = hb_start
    grp._v_attrs.z_scale              = 1./hb_dz
    grp._v_attrs.left_right_node      = 0.0
    grp._v_attrs.left_right_sharpness = 1./hb_dz
    grp._v_attrs.half_thickness       = half_thickness

    grp._v_attrs.use_curvature        = int(use_curvature)
    grp._v_attrs.curvature_radius     = curvature_radius
    grp._v_attrs.curvature_sign       = curvature_sign


    infer_group = t.get_node('/input/potential/infer_H_O')
    d_res_id    = infer_group.donors.residue[:]
    a_res_id    = infer_group.acceptors.residue[:]
    res_id      = np.concatenate([d_res_id, a_res_id])
    n_donor     = infer_group.donors   .id.shape[0]
    n_acceptor  = infer_group.acceptors.id.shape[0]
    n_hb        = n_donor + n_acceptor

    hb_uesd_id = []
    if len(hbond_exclude_residues) > 0 :
        for i,r in enumerate(res_id):
            if not (r in hbond_exclude_residues):
                hb_uesd_id.append(i)
        hb_uesd_id = np.array(hb_uesd_id)
    else:
        hb_uesd_id = np.arange(n_hb)

    hb_type = np.ones(n_hb, dtype='i')
    hb_type[:n_donor] *= 0

    create_array(grp, 'hb_index', hb_uesd_id)
    create_array(grp,  'hb_type', hb_type)
    create_array(grp,    'coeff', hb_energy)

def write_membrane_potential4(fasta_seq, membrane_potential_fpath, membrane_thickness, membrane_exclude_residues, hbond_exclude_residues, use_curvature=False, curvature_radius=200, curvature_sign=1):

    half_thickness = membrane_thickness * 0.5

    with tb.open_file(membrane_potential_fpath) as lib:
        resnames       = lib.root.names[:]
        resname_to_num = dict([(aa.decode('ASCII'),i) for i,aa in enumerate(resnames)])
        cb_energy      = lib.root.cb_energy[:]
        hb_energy      = lib.root.hb_energy[:]
        icb_energy     = lib.root.icb_energy[:]
        ihb_energy     = lib.root.ihb_energy[:]
        cb_z_min       = lib.root._v_attrs.cb_z_min
        cb_z_max       = lib.root._v_attrs.cb_z_max
        hb_z_min       = lib.root._v_attrs.hb_z_min
        hb_z_max       = lib.root._v_attrs.hb_z_max

        cb_n_node      = cb_energy.shape[2]
        hb_n_node      = hb_energy.shape[2]

        n_bl           = cb_energy.shape[1]
        burial_nodes   = lib.root.burial_nodes[:]
        assert (burial_nodes.size+1 == n_bl)

        cb_dz = (cb_z_max-cb_z_min)/(cb_n_node-2)
        hb_dz = (hb_z_max-hb_z_min)/(hb_n_node-2)

        cb_start = cb_z_min - half_thickness
        hb_start = hb_z_min - half_thickness


    curvature_center = Const3DCoord1(t, 'curvature_center', 'pos', 0, 1, 2, [[0., 0., curvature_radius*-1*curvature_sign ]])

    # memb potential for CB atoms
    grp = t.create_group(t.root.input.potential, 'cb_surf_membrane_potential')
    grp._v_attrs.arguments            = np.array([b'placement_fixed_point_only_CB', b'environment_coverage_sc', b'surface', bstring(curvature_center)])
    grp._v_attrs.integrator_level     = 0
    grp._v_attrs.z_start              = cb_start
    grp._v_attrs.z_scale              = 1./cb_dz
    grp._v_attrs.bl_range_sharpness   = 1.0
    grp._v_attrs.left_right_node      = 0.0
    grp._v_attrs.left_right_sharpness = 1./cb_dz
    grp._v_attrs.half_thickness       = half_thickness

    grp._v_attrs.use_curvature        = int(use_curvature)
    grp._v_attrs.curvature_radius     = curvature_radius
    grp._v_attrs.curvature_sign       = curvature_sign

    n_res = len(fasta_seq)
    used_res = []

    if len(membrane_exclude_residues) > 0:
        for r in range(n_res):
            if not (r in membrane_exclude_residues):
                used_res.append(r)
        used_res = np.array(used_res)
    else:
        used_res = np.arange(n_res)

    res_type = []
    for r in used_res:
        res_type.append(resname_to_num[fasta_seq[r]])
    res_type = np.array(res_type)

    env_group  = t.get_node('/input/potential/sigmoid_coupling_environment')
    weights    = env_group.weights[:20]

    create_array(grp,    'cb_index', used_res)
    create_array(grp,    'res_type', res_type)
    create_array(grp,     'weights', weights)
    create_array(grp,    'bl_nodes', burial_nodes)
    create_array(grp,       'coeff', cb_energy)
    create_array(grp, 'coeff_inner', icb_energy)

    # memb potential for HB atoms
    grp = t.create_group(t.root.input.potential, 'hb_surf_membrane_potential')
    grp._v_attrs.arguments            = np.array([b'protein_hbond', b'environment_coverage_sc', b'surface', bstring(curvature_center)])
    grp._v_attrs.integrator_level     = 0
    grp._v_attrs.z_start              = hb_start
    grp._v_attrs.z_scale              = 1./hb_dz
    grp._v_attrs.left_right_node      = 0.0
    grp._v_attrs.left_right_sharpness = 1./hb_dz
    grp._v_attrs.half_thickness       = half_thickness

    grp._v_attrs.use_curvature        = int(use_curvature)
    grp._v_attrs.curvature_radius     = curvature_radius
    grp._v_attrs.curvature_sign       = curvature_sign

    infer_group = t.get_node('/input/potential/infer_H_O')
    d_res_id    = infer_group.donors.residue[:]
    a_res_id    = infer_group.acceptors.residue[:]
    res_id      = np.concatenate([d_res_id, a_res_id])
    n_donor     = infer_group.donors   .id.shape[0]
    n_acceptor  = infer_group.acceptors.id.shape[0]
    n_hb        = n_donor + n_acceptor

    hb_uesd_id  = []
    hb_uesd_res = []
    if len(hbond_exclude_residues) > 0 :
        for i,r in enumerate(res_id):
            if not (r in hbond_exclude_residues):
                hb_uesd_id.append(i)
                hb_uesd_res.append(r)
        hb_uesd_id  = np.array(hb_uesd_id)
        hb_uesd_res = np.array(hb_uesd_res)
    else:
        hb_uesd_id = np.arange(n_hb)
        hb_uesd_res = res_id

    hb_type = np.ones(n_hb, dtype='i')
    hb_type[:n_donor] *= 0

    create_array(grp,    'hb_index', hb_uesd_id)
    create_array(grp,     'hb_type', hb_type)
    create_array(grp,      'hb_res', hb_uesd_res)
    create_array(grp,       'coeff', hb_energy)
    create_array(grp,     'weights', weights)
    create_array(grp, 'coeff_inner', ihb_energy)

def write_membrane_lateral_potential( parser, fasta_seq, lateral_potential_fpath):

    n_res = len(fasta_seq)
    resid = np.arange(n_res)

    grp = t.create_group(t.root.input.potential, 'membranelateral_potential')
    grp._v_attrs.arguments = np.array([b'placement_fixed_point_only_CB', b'environment_coverage_sc', b'surface'])

    #<----- ----- ----- ----- make energy splines ----- ----- ----- ----->#
    import scipy.interpolate
    def extrapolated_spline(x0, y0):
        spline = scipy.interpolate.InterpolatedUnivariateSpline(x0,y0)
        def f(x, spline=spline):
            return np.select(
                    [(x<x0[0]),              (x>x0[-1]),              np.ones_like(x,dtype='bool')],
                    [np.zeros_like(x)+y0[0], np.zeros_like(x)+y0[-1], spline(x)])
        return f

    fields = [ln.split() for ln in open(lateral_potential_fpath, 'r')]
    header_fields = 'z lateral_pressure'.split()
    if [x.lower() for x in fields[0]] == header_fields:
        len_header_fields = len(header_fields)
    else:
        parser.error('First line must be "%s"'%(" ".join(header_fields)))
    if not all(len(f)==len_header_fields for f in fields):
        parser.error('Invalid format')
    fields = fields[1:]
    n_point = len(fields)
    fields = np.array(fields)
    fields = fields.astype(float)

    min_z = np.min(fields[:,0])
    max_z = np.max(fields[:,0])
    n_knot = int(max_z-min_z)+1
    force_splines = extrapolated_spline(fields[:,0], fields[:,1])
    z_ = np.linspace(min_z, max_z, n_knot)
    f_ = force_splines(z_)
    
    import upside_engine as ue
    pp = ue.clamped_spline_solve(f_)

    #<----- ----- ----- ----- write to grp ----- ----- ----- ----->#
    cov_midpoint  = np.ones(20) + 2.5
    cov_sharpness = np.ones(20)

    if 'cb_surf_membrane_potential' in t.root.input.potential:
        used_res = t.root.input.potential.cb_surf_membrane_potential.cb_index[:]
        res_type = t.root.input.potential.cb_surf_membrane_potential.res_type[:]
        weights  = t.root.input.potential.cb_surf_membrane_potential.weights[:]
        half_thickness = t.root.input.potential.cb_surf_membrane_potential._v_attrs.half_thickness
    elif 'cb_membrane_potential' in t.root.input.potential:
        used_res = t.root.input.potential.cb_membrane_potential.cb_index[:]
        res_type = t.root.input.potential.cb_membrane_potential.res_type[:]
        weights  = t.root.input.potential.cb_membrane_potential.weights[:]
        half_thickness = t.root.input.potential.cb_membrane_potential._v_attrs.half_thickness

    create_array(grp,       'cb_index', used_res)
    create_array(grp,       'restypes', res_type)
    create_array(grp,        'weights', weights)
    create_array(grp,   'cov_midpoint', cov_midpoint)
    create_array(grp,  'cov_sharpness', cov_sharpness)
    create_array(grp,         'params', pp )

    grp._v_attrs.min_z = min_z
    grp._v_attrs.max_z = max_z
    grp._v_attrs.inv_dz = 1./(max_z-min_z)
    grp._v_attrs.n_knot = len(pp)
    grp._v_attrs.half_thickness = half_thickness
    grp._v_attrs.integrator_level     = 0

#---------------------------------------------------------------------------
#                      Utility functions for potentials
#---------------------------------------------------------------------------
def apply_param_scale(hb_scale=1., env_scale=1., rot_scale=1., memb_scale=1.):
    pot_group = t.root.input.potential
    if hb_scale != 1.:
        print ("scaling hb {}x".format(hb_scale))
        # See write_short_hbond() for layout of params.
        # WARNING: First three are supposed to be helix_hbond_energy, sheet_hbond_energy, turn_hbond_energy, 
        # but there seems to be an extra fourth energy in the actual param data not included in the (outdated?) layout  
        pot_group.hbond_energy.parameters[:4] *= hb_scale
        # ToDo: what about the hbond_coverage groups -> how are they coupled
        # to the energy? For env it's clear, but not so for hb (not sure if this still applies to Upside2)

    if env_scale != 1.:
        print ("scaling env {}x".format(env_scale))
        if "sigmoid_coupling_environment" in pot_group:
            pot_group.sigmoid_coupling_environment.scale[:] *= env_scale
        elif "nonlinear_coupling_environment" in pot_group:
            pot_group.nonlinear_coupling_environment.coeff[:] *= env_scale

        # Backbone environment
        if "bb_sigmoid_coupling_environment" in pot_group:
            pot_group.bb_sigmoid_coupling_environment._v_attrs.scale *= env_scale

    if rot_scale != 1.:
        print ("scaling rot {}x".format(rot_scale))
        # ToDo: scale placement energy? They are rama dep 1-body
        # Scaling only the radial parts
        n_knot_ang = 15
        pot_group.hbond_coverage.interaction_param[2*n_knot_ang:] *= rot_scale
        pot_group.hbond_coverage_hydrophobe.interaction_param[2*n_knot_ang:] *= rot_scale
        pot_group.rotamer.pair_interaction.interaction_param[2*n_knot_ang:] *= rot_scale

    if memb_scale != 1.:
        print ("scaling memb {}x".format(memb_scale))
        # Regular membrane potential is written by write_membrane_potential3()
        if "cb_membrane_potential" in pot_group:
            pot_group.cb_membrane_potential.coeff[:]  *= memb_scale
            pot_group.hb_membrane_potential.coeff[:] *= memb_scale
        elif "cb_surf_membrane_potential" in pot_group:
            pot_group.cb_surf_membrane_potential.coeff[:]  *= memb_scale
            pot_group.hb_surf_membrane_potential.coeff[:] *= memb_scale
    print()

#---------------------------------------------------------------------------
#                              main function 
#---------------------------------------------------------------------------


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Prepare input file',
            usage='use "%(prog)s --help" for more information')
    parser.add_argument('--fasta', required=True,
            help='[required] FASTA sequence file')
    parser.add_argument('--output', default='system.h5', required=True,
            help='path to output the created .h5 file (default system.h5)')
    parser.add_argument('--target-structure', default='',
            help='Add target .initial.npy structure for later analysis.  This information is written under '+
            '/target and is never read by Upside.  The /target group may be useful for later analysis.')
    parser.add_argument('--initial-structure', default='',
            help='Numpy pickled file for initial structure for the simulation.  ' +
            'If there are not enough structures for the number of replicas ' +
            'requested, structures will be recycled.  If not provided, a ' +
            'freely-jointed chain with good bond lengths and angles but bad dihedrals will be used ' +
            'instead.')
    parser.add_argument('--chain-break-from-file', default='',
            help='File with indices of chain first residues recorded during initial structure generation to automate creation and analysis of separate chains.')
    parser.add_argument('--jump-length-scale', type=float, default=5., help='Translational gaussian width in angstroms for Monte Carlo JumpSampler. Default: 5 angstroms')
    parser.add_argument('--jump-rotation-scale', type=float, default=30., help='Rotational gaussian width in degrees for Monte Carlo JumpSampler. Default: 30 degrees')
    parser.add_argument('--remove-pivot', action='store_true', help='Whether to remove the MC PivotSampler param group to isolate JumpSampler for testing')

    parser.add_argument('--bond-stiffness', default=48., type=float,
            help='Bond spring constant in units of energy/A^2 (default 48)')
    parser.add_argument('--angle-stiffness', default=175., type=float,
            help='Angle spring constant in units of 1/dot_product (default 175)')
    parser.add_argument('--omega-stiffness', default=30., type=float,
            help='Omega angle spring constant in units of 1/dot_product (default 30)')
    parser.add_argument('--debugging-only-disable-basic-springs', default=False, action='store_true',
            help='Disable basic springs (like bond distance and angle).  Do not use this.')

    parser.add_argument('--intensive-memory', default=False, action='store_true',
            help='Avoid sparse matrices to significantly save memory consumption. This option is suitable for large ' +
            'simulation systems, such as > 1000 amino acids.')
    parser.add_argument('--no-intensive-memory', default=False, action='store_true',
            help='Turn off the intensive memory mode.')

    parser.add_argument('--rama-library', default='',
            help='smooth Rama probability library')
    parser.add_argument('--rama-library-combining-rule', default='mixture',
            help='How to combine left and right coil distributions in Rama library '+
            '(mixture or product).  Default is mixture.')
    parser.add_argument('--trans-cis', default='', 
            help='Table of the trans-cis state for PRO.  Each line must contain 4 fields and the header line must be '+
            '"residue init_state energy k width". '+
            '"residue" is used to specify the residue id of PRO involved in trans-cis transformation (starting from 0). ' +
            '"init_state" is the initial state (0: trans, 1: cis). ' +
            '"energy" is the energy barrier from trans state to the transition state. ' +
            '"k" is used to set the  energy barrier from cis state to the transition state: k*energy. ' +
            '"width" controls the 1/4 width (degree) of the energy barrier (such as 45 degree).')
    parser.add_argument('--rama-sheet-library', default=None,
            help='smooth Rama probability library for sheet structures')
    parser.add_argument('--secstr-bias', default='',
            help='Bias file for secondary structure.  First line of the file must be "residue secstr energy".  '+
            'secstr must be one of "helix" or "sheet".  Bias is implemented by a simple Rama bias, hence coil bias '+
            'is not implemented.')
    parser.add_argument('--rama-sheet-mixing-energy', default='', 
            help='reference energy for sheets when mixing with coil library.  More negative numbers mean more '+
            'sheet content in the final structure.  Default is no sheet mixing.')
    parser.add_argument('--rama-param-deriv', default=False, action='store_true',
            help='generate the param deriv for rama potential')
    parser.add_argument('--reference-state-rama', default='',
            help='Do not use this unless you know what you are doing.')

    parser.add_argument('--no-backbone', dest='backbone', default=True, action='store_false',
            help='do not use rigid nonbonded for backbone N, CA, C, and CB')

    parser.add_argument('--hbond-energy', default=0.,
            help='energy for forming a protein-protein hydrogen bond.  Default is no HBond energy.')
    parser.add_argument('--hbond-exclude-residues', default=[], type=parse_segments,
            help='Residues to have neither hydrogen bond donors or acceptors')
    parser.add_argument('--loose-hbond-criteria', default=False, action='store_true',
            help='Use far more permissive angles and distances to judge HBonding.  Do not use for simulation. '+
            'This is only useful for static backbone training when crystal or NMR structures have poor '+
            'hbond geometry.')

    parser.add_argument('--rotamer-placement', default=None,
            help='rotameric sidechain library')
    parser.add_argument('--dynamic-rotamer-placement', default=False, action='store_true',
            help='Use dynamic rotamer placement (not recommended)')
    parser.add_argument('--dynamic-rotamer-1body', default=False, action='store_true',
            help='Use dynamic rotamer 1body')
    parser.add_argument('--fix-rotamer', default='',
            help='Table of fixed rotamers for specific sidechains.  A header line must be present and the first '+
            'three columns of that header must be '+
            '"residue restype rotamer", but there can be additional, ignored columns.  The restype must '+
            'match the corresponding restype in the FASTA file (intended to prevent errors).  It is permissible '+
            'to fix only a subset of rotamers.  The value of rotamer must be an integer, corresponding to the '+
            'numbering in the --rotamer-placement file.  Such a file can be created with PDB_to_initial_structure '+
            '--output-chi1.')
    parser.add_argument('--rotamer-interaction', default=None,
            help='rotamer sidechain pair interaction parameters')
    parser.add_argument('--rotamer-solve-damping', default=0.4, type=float,
            help='damping factor to use for solving sidechain placement problem')
    parser.add_argument('--rotamer-exclude-residues', default=[], type=parse_segments,
            help='Residues that do not have rotamer sidechain beads nor participate in their pair interaction. Currently not properly implemented')

    parser.add_argument('--environment-potential', default='',
            help='Path to many-body environment potential')
    parser.add_argument('--environment-potential-type', default=1, type=int,
            help='The type 0 is spline-like nonlinear coupling function, 1 is the sigmoid-like coupling function')
    parser.add_argument('--environment-weights-number', default=20, type=int,
            help='Only three values are allowed: 1, 20, 40. The value 1 means that the weights of the side chain ' + 
            'beads of all amino acids are the same (1.0); The value 20 means that the weights depends on the amino ' +
            'acid type of side chain beads; The value 400 means that the weights depends on the residue type type of ' + 
            'side chain beads or CB atoms belongs to.')
    parser.add_argument('--vector-CA-CO', default=False, action='store_true', help='use vector CA-CO for the hemisphere orientation to be compatible with FF2.0 ')
    parser.add_argument('--env-exclude-residues', default=[], type=parse_segments,
            help='Residues that do not participate as the CB center for the environment interaction. Note that residues ' +
            'in --rotamer-exclude-residues do not contribute beads for the interaction. Currently not properly implemented')

    parser.add_argument('--bb-environment-potential', default='',
            help='Path to many-body environment potential for backbone')
    parser.add_argument('--use-heavy-atom-coverage', default=False, action='store_true',
            help='include heavy atoms on backbone when count the coverage for H/O')

    parser.add_argument('--surface', default=False, action='store_true', help='surface')
    parser.add_argument('--surface-method', default=1, type=int, help='surface finding method. 0: fast method; 1: lipid diffusion algorithm.')
    parser.add_argument('--surface-included-residues', default=[], type=parse_segments,
            help='List of residues in the protein.  The residue list should be of a form like ' +
            '--restraint-group 10-13,17,19-21 and that list would specify all the atoms in '+
            'residues 10,11,12,13,17,19,20,21. '+
            'Each atom in the specified residues will be randomly connected to atoms in other residues by ' +
            'springs with equilibrium distance given by the distance of the atoms in the initial structure.  ' +
            'Multiple restraint groups may be specified by giving the --restraint-group flag multiple times '
            'with different residue lists.  The strength of the restraint is given by --restraint-spring-constant')
    parser.add_argument('--zbuffer', default=None, type=float, help='buffer region to determine the residue ' +
            'list which is subject to the membrane force field.')

    parser.add_argument('--membrane-thickness', default=None, type=float,
            help='Thickness of the membrane in angstroms for use with --membrane-potential.')
    parser.add_argument('--channel-membrane-potential', default='',
            help='Parameter file (.h5 format) for membrane potential with channel information. User must also supply --membrane-thickness.')
    parser.add_argument('--membrane-potential', default='',
            help='Parameter file (.h5 format) for membrane potential. User must also supply --membrane-thickness.')

    parser.add_argument('--use-curvature', default=False, action='store_true', help='model the curvature layers of lipid.')
    parser.add_argument('--curvature-radius', default=1000., type=float,
            help='Curvature radius of the membrane in angstroms for use with --membrane-potential and --use-curvature.')
    parser.add_argument('--curvature-sign', default=1., type=int,
            help='Curvature sign of the membrane (1 for positive curvature,  -1 for negative curvature)')

    parser.add_argument('--membrane-exclude-residues', default=[], type=parse_segments,
            help='Residues that do not participate in the --membrane-potential (same format as --restraint-group).' +
                 'User must also supply --membrane-potential.')
    parser.add_argument('--membrane-exposed-criterion', default=None,
            help='Parameter file including the sigmoid-like exposed criterion (.h5 format) for membrane potential.')

    parser.add_argument('--membrane-lateral-potential', default='',
            help='Table of lateral pressure curve.  Each line must contain 2 fields and the header line must be '+
            '"z lateral_pressure". Next is followed by the z position and the corresponding pressure value in bar. '+
            'The largest z and the smallest z will be the boundaries of the lateral pressure.')

    parser.add_argument('--hb-scale', default=1., type=float, help='Scaling for hbond potential')
    parser.add_argument('--env-scale', default=1., type=float, help='Scaling for environment potentials')
    parser.add_argument('--rot-scale', default=1., type=float, help='Scaling for rotamer pair potential')
    parser.add_argument('--memb-scale', default=1., type=float, help='Scaling for membrane potential')

    args = parser.parse_args()

    if args.rotamer_exclude_residues or args.env_exclude_residues:
        parser.error('--rotamer-exclude-residues and --env-exclude-residues are not properly implemented yet')

    require_affine = False
    require_rama = False
    require_backbone_point = False
    require_weighted_pos = False

    fasta_seq_with_cpr = read_fasta(open(args.fasta))
    fasta_seq = np.array([(x if x != 'CPR' else 'PRO') for x in fasta_seq_with_cpr])  # most potentials don't care about CPR

    global n_atom, t, potential, n_chains, chain_first_residue, chain_starts, rl_chains, use_intensive_memory

    n_res = len(fasta_seq)
    n_atom = 3*n_res

    use_intensive_memory = False
    if n_res > 1000:
        use_intensive_memory = True

    t = tb.open_file(args.output,'w')

    input = t.create_group(t.root, 'input')
    create_array(input, 'sequence', obj=fasta_seq_with_cpr)


    # the initial structure:
    #---------------------------------------------

    if args.initial_structure:
        init_pos = np.load(args.initial_structure)
        assert init_pos.shape == (n_atom, 3)

    if args.target_structure:
        def f():
            # little function closure to protect the namespace from ever seeing the target structure
            target_pos = np.load(args.target_structure)
            assert target_pos.shape == (n_atom, 3)
            g_target = t.create_group(t.root, 'target')
            target_pos = target_pos[:,:,None]
            t.create_array(t.root.target, 'pos', obj=target_pos[:,:,0])
        f()

    pos = np.zeros((n_atom, 3, 1), dtype='f4')
    if args.initial_structure:
        pos[:,:,0] = init_pos
    else:
        pos[:,:,0] = random_initial_config(len(fasta_seq))
    create_array(input, 'pos', obj=pos)

    if args.intensive_memory:
        use_intensive_memory = True
    if args.no_intensive_memory:
        use_intensive_memory = False


    # the chains:
    #---------------------------------------------

    n_chains = 1
    if args.chain_break_from_file:
        try:
            with open(args.chain_break_from_file) as infile:
                chain_dat = list(infile)
            # chain_first_residue = np.loadtxt(args.chain_break_from_file, ndmin=1, dtype='int32')
        except IOError:
            chain_first_residue = np.array([], dtype='int32')
            n_chains = 1
        else:
            if len(chain_dat) > 2:
                has_rl_info = True
            else:
                has_rl_info = False
            chain_first_residue = np.array(chain_dat[0].split(), dtype='int32')
            chain_counts = np.array(chain_dat[1].split(), dtype='int32')
            print ("chain_first_residue:", chain_first_residue)
            print ("chain_counts:", chain_counts)
            n_chains = chain_first_residue.size+1
            if has_rl_info:
                rl_chains_actual = np.array(chain_dat[2].split(), dtype='int32')
                rl_chains = np.array(chain_dat[3].split(), dtype='int32')
                print ("rl_chains_actual:", rl_chains_actual)
                print ("rl_chains:", rl_chains)
            else:
                rl_chains = None

        print ()
        print ("n_chains:", n_chains)

        if chain_first_residue.size:
            break_grp = t.create_group("/input","chain_break","Indicates that multi-chain simulation and removal of bonded potential terms accross chains requested")
            t.create_array(break_grp, "chain_first_residue", chain_first_residue, "Contains array of chain first residues, apart from residue 0")
            t.create_array(break_grp, "chain_counts", chain_counts, "Counts of (broken) chains for each actual chain")
            if has_rl_info:
                t.create_array(break_grp, "rl_chains_actual", rl_chains_actual, "Numbers of actual receptor and ligand chains")
                t.create_array(break_grp, "rl_chains", rl_chains, "Numbers of receptor and ligand chains")    

            required_hbond_exclude_res = [i+j for i in chain_first_residue for j in [-1,0]]
            if args.hbond_exclude_residues:
                args.hbond_exclude_residues = np.unique(np.append(args.hbond_exclude_residues, required_hbond_exclude_res))
            else:
                args.hbond_exclude_residues = np.array(required_hbond_exclude_res)

            print ()
            print ("hbond_exclude_residues")
            print (args.hbond_exclude_residues)

        chain_starts = np.array(chain_first_residue)*3
        chain_starts = np.append([0], chain_starts)

    # record the commands and options
    #---------------------------------------------

    args_group = t.create_group(input, 'args')
    for k,v in sorted(vars(args).items()):
        args_group._v_attrs[k] = v
    args_group._v_attrs['invocation'] = ' '.join(sys.argv[:])


    # potential 
    #---------------------------------------------

    potential = t.create_group(input,  'potential')

    # bonded terms

    if not args.debugging_only_disable_basic_springs:
        write_dist_spring(args)
        write_angle_spring(args)
        if args.trans_cis:
            write_omega_spring2(parser, args, fasta_seq_with_cpr, args.trans_cis)
        else:
            write_omega_spring1(args, fasta_seq_with_cpr)

    if args.rama_library:
        require_rama = True
        if args.trans_cis:
            write_rama_map_pot2(parser, fasta_seq_with_cpr, args.rama_library, args.trans_cis, 
                                args.rama_sheet_mixing_energy, args.rama_library_combining_rule)
        else:
            write_rama_map_pot(fasta_seq_with_cpr, args.rama_library, args.rama_sheet_mixing_energy,
                args.secstr_bias, args.rama_library_combining_rule, args.rama_param_deriv)
    # elif args.torus_dbn_library:
    #     require_rama = True
    #     write_torus_dbn(fasta_seq_with_cpr, args.torus_dbn_library)
    else:
        eprint('WARNING: running without any Rama potential !!!')


    # hack to fix reference state issues for Rama potential
    if args.reference_state_rama:
        # define correction
        ref_state_cor =  np.log(cPickle.load(open(args.reference_state_rama, 'rb'), encoding='latin1'))
        ref_state_cor -= ref_state_cor.mean()

        grp = t.create_group(potential, 'rama_map_pot_ref')
        grp._v_attrs.arguments = np.array([b'rama_coord'])
        grp._v_attrs.log_pot = 0

        create_array(grp, 'residue_id',   obj=np.arange(len(fasta_seq)))
        create_array(grp, 'rama_map_id',  obj=np.zeros(len(fasta_seq), dtype='i4'))
        create_array(grp, 'rama_pot',     obj=ref_state_cor[None])

    # nonbonded terms

    if args.backbone:
        require_affine = True
        write_backbone_pair(fasta_seq)

    if args.hbond_energy:
        write_infer_H_O  (fasta_seq, args.hbond_exclude_residues)
        write_count_hbond(fasta_seq, args.loose_hbond_criteria )
        write_short_hbond(fasta_seq, args.hbond_energy)#, args.hbond_energy_patterns)

    sc_node_name = ''
    if args.rotamer_placement:
        require_rama = True
        require_affine = True
        sc_node_name, pl_node_name = write_rotamer_placement(
                fasta_seq, args.rotamer_placement,
                args.dynamic_rotamer_placement, args.dynamic_rotamer_1body,
                args.fix_rotamer, args.rotamer_exclude_residues)

    if args.hbond_energy and sc_node_name:
        write_rotamer_backbone(fasta_seq, args.rotamer_interaction, sc_node_name)

    if args.rotamer_interaction:
        # must be after write_count_hbond if hbond_coverage is used
        write_rotamer(fasta_seq, args.rotamer_interaction, args.rotamer_solve_damping, sc_node_name, pl_node_name)


    # multibody terms

    if args.environment_potential:
        if args.rotamer_placement is None:
            parser.error('--rotamer-placement is required, based on other options.')
        require_weighted_pos = True     
        write_environment(fasta_seq, args.environment_potential, sc_node_name, args.environment_potential_type, args.environment_weights_number, args.vector_CA_CO, False, args.env_exclude_residues, args.rotamer_exclude_residues)

    if args.bb_environment_potential:
        if args.environment_potential is None:
            parser.error('--environment-potential is required, based on other options.')
        require_weighted_pos = True     
        write_bb_environment(fasta_seq, args.environment_potential, sc_node_name, args.bb_environment_potential, args.use_heavy_atom_coverage)


    # membrance potential

    if args.surface:
        require_backbone_point = True
        if args.surface_included_residues:
            write_surface_coord(fasta_seq, args.surface_method, thickness=args.membrane_thickness, included_list=args.surface_included_residues )
        else:
            write_surface_coord(fasta_seq, args.surface_method, thickness=args.membrane_thickness)

    if args.channel_membrane_potential:
        if args.membrane_thickness is None:
            parser.error('--channel-membrane-potential requires --membrane-thickness')
        if args.surface is None:
            parser.error('--channel-membrane-potential requires --surface')
        require_backbone_point = True
        write_membrane_potential4(fasta_seq,
                                 args.channel_membrane_potential,
                                 args.membrane_thickness,
                                 args.membrane_exclude_residues, 
                                 args.hbond_exclude_residues,
                                 args.use_curvature,
                                 args.curvature_radius,
                                 args.curvature_sign)

    if args.membrane_potential:
        if args.membrane_thickness is None:
            parser.error('--membrane-potential requires --membrane-thickness')
        require_backbone_point = True
        write_membrane_potential3(fasta_seq,
                                 args.membrane_potential,
                                 args.membrane_thickness,
                                 args.membrane_exclude_residues, 
                                 args.hbond_exclude_residues,
                                 args.use_curvature,
                                 args.curvature_radius,
                                 args.curvature_sign)

    if args.membrane_lateral_potential:
        if not args.surface:
            parser.error('--membrane-lateral-potential requires --surface')
        if args.membrane_thickness is None:
            parser.error('--membrane-lateral-potential requires --membrane-thickness')
        if args.membrane_potential=='' and args.channel_membrane_potential=='':
            parser.error('--membrane-lateral-potential requires --membrane-potential or --channel-membrane-potential')
        write_membrane_lateral_potential( parser, fasta_seq, args.membrane_lateral_potential)

    if require_backbone_point:
        require_affine = True
        write_CB(fasta_seq)

    if require_weighted_pos:
        write_weighted_pos(sc_node_name, pl_node_name)

    if require_rama:
        #write_rama_coord()
        write_rama_coord2()

    if require_affine:
        write_affine_alignment(len(fasta_seq))


    # if we have the necessary information, write pivot_sampler
    if require_rama and 'rama_map_pot' in potential:

        if n_chains > 1:
            # Need to add one atom past the last atom so that the last chain is processed
            chain_starts_plus = np.append(chain_starts, n_atom)

            if rl_chains is None:
                # range covers each individual chain
                jump_atom_range = np.array([[chain_starts_plus[i], chain_starts_plus[i+1]] for i in range(n_chains)], dtype='int32')
            else:
                # Add ranges of all chains in receptor and ligand for collective jumps
                jump_atom_range = np.array([[chain_starts_plus[0], chain_starts_plus[rl_chains[0]]]])
                jump_atom_range = np.append(jump_atom_range, [[chain_starts_plus[rl_chains[0]], chain_starts_plus[-1]]], axis=0)

            jump_sigma_trans = np.array([args.jump_length_scale]*len(jump_atom_range), dtype='float32')
            jump_sigma_rot   = np.array([args.jump_rotation_scale*np.pi/180.]*len(jump_atom_range), dtype='float32') # Converts to radians

            print ()
            print ("jump atom_range:\n{}\nsigma_trans:\n{}\nsigma_rot:\n{}\n".format(jump_atom_range, jump_sigma_trans, jump_sigma_rot))

            jump_grp = t.create_group("/input","jump_moves","JumpSampler Params")
            t.create_array(jump_grp, "atom_range",  jump_atom_range,  "First, last atom num demarking each chain")
            t.create_array(jump_grp, "sigma_trans", jump_sigma_trans, "Translational gaussian width")
            t.create_array(jump_grp, "sigma_rot",   jump_sigma_rot,   "Rotational gaussian width")
        else:
            grp = t.create_group(input, 'pivot_moves')
            pivot_atom = potential.rama_coord.id[:]
            non_terminal_residue = np.array([not(np.int64(-1).astype(pivot_atom.dtype) in tuple(x))
                for x in pivot_atom])

            create_array(grp, 'proposal_pot',  potential.rama_map_pot.rama_pot[:])
            create_array(grp, 'pivot_atom',    pivot_atom[non_terminal_residue])
            create_array(grp, 'pivot_restype', potential.rama_map_pot.rama_map_id_all[:][non_terminal_residue])
            create_array(grp, 'pivot_range',   np.column_stack((grp.pivot_atom[:,4]+1,np.zeros(sum(non_terminal_residue),'i')+n_atom)))

    # Scale potential if requested
    apply_param_scale(args.hb_scale, args.env_scale, args.rot_scale, args.memb_scale)

    t.close()

if __name__ == '__main__':
    main()
