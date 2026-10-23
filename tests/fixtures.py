#!/usr/bin/env python
# Build the audit fixtures: ff_2.1 configs of a few example proteins, with and without the membrane, then one
# config per option group of upside_config and advanced_config (a protein with one option, or with a few
# compatible ones).  The option tables are written beside the configs.
# usage: fixtures.py OUTPUT_DIR
import sys, os
import subprocess as sp
import numpy as np
import tables as tb
from scipy.spatial.transform import Rotation

upside_path = os.environ['UPSIDE_HOME']
sys.path.insert(0, upside_path + '/py')
import run_upside as ru

out_dir = sys.argv[1]
param_dir = upside_path + '/parameters/'
ff_dir = param_dir + 'ff_2.1/'

chig = 'example/01.GettingStarted/pdb/chig.pdb'
dfn  = 'example/01.GettingStarted/pdb/1dfn.pdb'
ubq  = 'example/07.MoreRestraints/pdb/1UBQ.pdb'
qhj  = 'example/08.MembraneSimulation/pdb/1qhj.pdb'
rkl  = 'example/08.MembraneSimulation/pdb/1rkl.pdb'

def table(name, text):
    path = os.path.join(out_dir, name)
    with open(path, 'w') as f: f.write(text)
    return path

def hbond_with_classes(name):
    # ff_2.1's H-bond energies with offsets for two residue classes; no shipped file has classes
    path = os.path.join(out_dir, name)
    with tb.open_file(ff_dir + 'hbond.h5') as t:
        params = t.root.parameter[:]
    with tb.open_file(path, 'w') as t:
        t.create_array(t.root, 'parameter', np.concatenate([params, [0.3, -0.2, 0.1, -0.1, 0.2, -0.3]]))
        t.create_array(t.root, 'class_restype', np.array([b'GLY', b'THR']))
    return path

def turn_omega(name, residue, angle):
    # the chain from N(residue) on, turned about the C(residue-1)-N(residue) bond: only omega(residue) changes;
    # no example structure has a cis proline
    x = np.load(os.path.join(out_dir, name + '.initial.npy'))
    c, n = x[3*residue-1], x[3*residue]
    x[3*residue:] = Rotation.from_rotvec(np.deg2rad(angle)*(n-c)/np.linalg.norm(n-c)).apply(x[3*residue:]-n) + n
    path = os.path.join(out_dir, '%s.omega%d.npy' % (name, residue))
    np.save(path, x)
    return path

def build(name, pdb, cli=[], adv={}, **options):
    base = os.path.join(out_dir, name)
    sp.check_call([sys.executable, upside_path + '/py/PDB_to_initial_structure.py',
                   os.path.join(upside_path, pdb), base, '--record-chain-breaks', '--disable-recentering'])
    kwargs = dict(rama_library=param_dir + 'common/rama.dat',
                  rama_sheet_mix_energy=ff_dir + 'sheet',
                  reference_state_rama=param_dir + 'common/rama_reference.pkl',
                  hbond_energy=ff_dir + 'hbond.h5',
                  rotamer_placement=ff_dir + 'sidechain.h5',
                  dynamic_rotamer_1body=True,
                  rotamer_interaction=ff_dir + 'sidechain.h5',
                  environment_potential=ff_dir + 'environment.h5',
                  bb_environment_potential=ff_dir + 'bb_env.dat',
                  chain_break_from_file=base + '.chain_breaks',
                  initial_structure=base + '.initial.npy')
    kwargs.update(options)
    command = ru.upside_config(base + '.fasta', base + '.up', **kwargs).split('\n')[0].split()
    if cli:  # options of upside_config.py that run_upside.upside_config does not pass on
        sp.check_call(command + cli)
    if adv:
        ru.advanced_config(base + '.up', **adv)
    print(base + '.up')

os.makedirs(out_dir, exist_ok=True)
membrane = dict(membrane_potential=ff_dir + 'membrane.h5', membrane_thickness=24.8)
centre_1dfn = table('1dfn.centre', 'resid1 resid2 axis_x axis_y axis_z offset_x offset_y offset_z\n'
                                   '10 40 1 0 0 0 0 0\n')

build('chig', chig)
build('1dfn', dfn)
build('1ubq', ubq)
build('1qhj_memb', qhj, membrane_potential=ff_dir + 'membrane.h5', membrane_thickness=31.8)

# upside_config options
build('chig_springs', chig, bond_stiffness=40., angle_stiffness=150., omega_stiffness=25., backbone=False)
build('chig_rama', chig, rama_param_deriv=True,
      cli=['--rama-library-combining-rule=product',
           '--secstr-bias=' + table('chig_rama.secstr', 'residue secstr energy\n1 sheet -0.5\n4 helix -0.5\n8 sheet -0.5\n')])
build('chig_transcis', chig, trans_cis=table('chig_transcis.pro', 'residue energy k width\n3 3 0.7 45.0\n'))
build('chig_transcis_cis', chig, trans_cis=os.path.join(out_dir, 'chig_transcis.pro'),
      initial_structure=turn_omega('chig', 3, 135.))
build('chig_hbond', chig, hbond_energy=hbond_with_classes('chig_hbond.hbond.h5'), hbond_exclude_residues='4',
      hb_scale=0.9, cli=['--loose-hbond-criteria'])
build('1ubq_rotamer', ubq, fix_rotamer=os.path.join(out_dir, '1ubq_rotamer.chi'), cli=['--rotamer-solve-damping=0.3'])
build('chig_env0', chig, environment_potential_type=0, vector_CA_CO=True)
build('chig_env0_w400', chig, environment_potential_type=0, vector_CA_CO=True, environment_weights_number=400)
build('chig_envw400', chig, environment_weights_number=400)
build('chig_hbbb', chig, use_heavy_atom_coverage=True, intensive_memory=True)
build('1rkl_curv', rkl, use_curvature=True, curvature_radius=120., curvature_sign=1, **membrane)
build('1rkl_scales', rkl, hb_scale=0.8, env_scale=1.2, rot_scale=0.9, memb_scale=1.1, **membrane)
build('1rkl_surface', rkl, surface=True, membrane_exclude_residues='0-2,28-30', **membrane)
build('1rkl_channel', rkl, surface=True, channel_membrane_potential=ff_dir + 'membrane.h5', membrane_thickness=24.8)
build('1rkl_channel_m0', rkl, surface=True, surface_method=0, surface_included_residues='4-26', zbuffer=2.,
      channel_membrane_potential=ff_dir + 'membrane.h5', membrane_thickness=24.8,
      membrane_exclude_residues='0-2', hbond_exclude_residues='29-30')

# advanced_config options
build('chig_restraints_xyz', chig, adv=dict(
      fixed_wall=table('chig_xyz.wall', 'residue radius spring_const wall_type x0 y0 z0\n'
                                        '0 3.0 4.0 1 -4.0 0.0 1.0\n9 2.0 4.0 0 -1.0 -3.0 3.0\n'),
      pair_wall=table('chig_xyz.pair_wall', 'residue1 residue2 radius spring_const wall_type xyz\n0 9 5.0 4.0 1\n'),
      fixed_spring=table('chig_xyz.spring', 'residue radius spring_const x0 y0 z0\n2-3 1.0 2.0 0.0 2.0 1.0\n'),
      pair_spring=table('chig_xyz.pair_spring', 'residue1 residue2 radius spring_const xyz\n1 8 5.0 2.0\n'),
      nail=table('chig_xyz.nail', 'residue spring_const xyz\n4-6 1.0\n')))
build('chig_restraints_1d2d', chig, adv=dict(
      fixed_wall=table('chig_1d2d.wall', 'residue radius spring_const wall_type x0\n0-1 2.0 4.0 1 -3.0\n5 1.0 4.0 0 4.0\n'),
      pair_wall=table('chig_1d2d.pair_wall', 'residue1 residue2 radius spring_const wall_type yz\n0 9 3.0 4.0 1\n'),
      fixed_spring=table('chig_1d2d.spring', 'residue radius spring_const x0 z0\n3-4 1.0 2.0 2.0 0.5\n'),
      pair_spring=table('chig_1d2d.pair_spring', 'residue1 residue2 radius spring_const z\n1 8 1.0 3.0\n'),
      nail=table('chig_1d2d.nail', 'residue spring_const xy\n2,6 1.0\n')))
build('chig_groups', chig, adv=dict(
      group=table('chig_groups.group', 'residues\n0-4\n5-9\n'),
      group_dist_spring=table('chig_groups.dist', 'group1 group2 radius spring_const xyz\n0 1 4.0 2.0\n'),
      group_nail_spring=table('chig_groups.nail', 'group spring_const y\n0 1.0\n'),
      group_const_spring=table('chig_groups.const', 'group radius spring_const x0 y0 z0\n1 1.0 2.0 0.0 0.0 0.0\n')))
build('chig_restraint', chig, adv=dict(restraint_groups=['0-4', '3-9'], restraint_spring_constant=2., cavity_radius=6.))
build('chig_contact', chig, adv=dict(
      contact_energies=table('chig_contact.contact', 'residue1 residue2 energy distance transition_width\n'
                                                     '1 8 -2.0 6.0 2.0\n2 7 -1.0 5.5 1.5\n'),
      tension=table('chig_contact.tension', 'residue tension_x tension_y tension_z\n0 0.1 0 0\n9 -0.1 0 0.05\n')))
build('chig_afm', chig, adv=dict(ask_before_using_AFM=table('chig_afm.afm',
      'residue spring_const tip_pos_x tip_pos_y tip_pos_z pulling_vel_x pulling_vel_y pulling_vel_z\n'
      '0 0.05 -9.0 -1.0 3.0 0.0 0.0 0.0\n9 0.05 -1.0 -6.0 6.0 0.0 0.0 0.001\n')))
build('chig_afm_auto', chig, adv=dict(ask_before_using_AFM=table('chig_afm_auto.afm',
      'residue spring_const pulling_vel_x pulling_vel_y pulling_vel_z\n0 0.05 0.0 0.0 0.0\n9 0.05 0.0 0.0 0.1\n'),
      AFM_time_step=0.009))
build('1dfn_unbound', dfn, spatial_transform_from_table=centre_1dfn, cli=['--jump-length-scale=3', '--jump-rotation-scale=20'],
      adv=dict(heuristic_cavity_radius=1.5, make_unbound=True, apply_restraint_group_to_each_chain=True))
build('1dfn_cavity', dfn, spatial_transform_from_table=centre_1dfn,
      adv=dict(cavity_radius_from_config=os.path.join(out_dir, '1dfn_unbound.up')))
