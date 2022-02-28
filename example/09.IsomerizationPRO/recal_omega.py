import os, sys
import numpy as np
import tables as tb
import mdtraj as md

upside_path = os.environ['UPSIDE_HOME']
upside_utils_dir = os.path.expanduser(upside_path+"/py")
sys.path.insert(0, upside_utils_dir)

import mdtraj_upside as mu
import upside_engine as ue

up_file   = sys.argv[1]
traj_file = sys.argv[2]

traj          = mu.load_upside_traj(traj_file)
pos           = mu.extract_bb_pos_angstroms(traj)[:]
engine        = ue.Upside(up_file)
n_frame       = traj.n_frames

for i in range(n_frame):
    engine.energy(pos[i])
    omega = engine.get_output('Dihedral_OmegaTransCis')[:,0]
    lamba = engine.get_output('SigmoidCoord_trans1')[:,0]
    lambb = engine.get_output('SigmoidCoord_trans2')[:,0]
    lamb1 = engine.get_output('Add_lambda_trans')[:,0]
    lamb2 = engine.get_output('Multiply_lambda_cis')[:,0]
    print (omega[0], lamb1[0], lamb2[0], lamba[0], lambb[0])
