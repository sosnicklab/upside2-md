import os, sys
import numpy as np
import tables as tb
import mdtraj as md

upside_path = os.environ['UPSIDE_HOME']
upside_utils_dir = os.path.expanduser(upside_path+"/py")
sys.path.insert(0, upside_utils_dir)

import mdtraj_upside as mu
import upside_engine as ue

pot = ["SigmoidEnergy_cis0", "SigmoidEnergy_cis1", "SigmoidEnergy_trans0", "SigmoidEnergy_trans1", "RamaMap2_cis", "RamaMap2_trans"]

up_file   = sys.argv[1]
traj_file = sys.argv[2]

traj          = mu.load_upside_traj(traj_file)
pos           = mu.extract_bb_pos_angstroms(traj)[:]
engine        = ue.Upside(up_file)
n_frame       = traj.n_frames
n_energy_term = len(pot)

energies = np.zeros([n_energy_term, n_frame])

for i in range(n_frame):
    engine.energy(pos[i])
    for j in range(n_energy_term):
        energies[j, i] = engine.get_output(pot[j])[0,0]

for i in range(n_frame):
    print (energies[0,i], energies[1,i], energies[2,i], energies[3,i], energies[4,i], energies[5,i])
