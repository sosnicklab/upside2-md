import sys, os, time, shutil
from math import sqrt
import subprocess as sp
import numpy as np

upside_path = '/home/pengxd/upside-ff2.0v/'
upside_utils_dir = os.path.expanduser(upside_path+"/py")
sys.path.insert(0, upside_utils_dir)
import run_upside as ru

## General Settings and Path
pdb_id = sys.argv[1]
#thickness = sys.argv[2]
#is_native = bool(int(sys.argv[2]))
is_native = False
is_native = True

# Parameters
param_dir_base = os.path.expanduser(upside_path+"/parameters/")
param_dir_common = param_dir_base + "common/"
param_dir_ff1 = param_dir_base + "ff_2.1/"

## Generate Upside readable initial structure (and fasta) from PDB 
print ("Initial structure gen...")
cmd = (
       "python {0}/PDB_to_initial_structure.py "
       "{1}.pdb "
       "{1} "
       "--record-chain-breaks "
       "--disable-recentering "
       "--allow-unexpected-chain-breaks "
      ).format(upside_utils_dir, pdb_id )
print (cmd)
sp.check_output(cmd.split())

### Configure
#print ("Configuring...")
#fasta = "{}.fasta".format(pdb_id)
#kwargs = dict(
#               rama_pot = param_dir_common + "rama.dat",
#               sheet_mix_energy = param_dir_ff1 + 'sheet',
#               hbond = param_dir_ff1 + "hbond.h5",
#               reference_rama = param_dir_common + "rama_reference.pkl",
#               placement = param_dir_ff1 + "sidechain.h5",
#               dynamic_1body = True,
#               rotamer_interaction_param = param_dir_ff1 + "sidechain.h5",
#               environment = param_dir_ff1 + "environment.h5",
#               bb_environment = param_dir_ff1 + "bb_env.dat",
#               environment_type = 1,
#               chain_break_from_file = "{}.chain_breaks".format(pdb_id),
#               #membrane_potential = "membrane_new.h5",
#               #membrane_thickness = float(thickness),
#               surface = True,
#               )
#
#if is_native:
#    kwargs['init'] =  "{}.initial.pkl".format(pdb_id)
#
#config_base = "{}.up".format(pdb_id)
#config_stdout = ru.upside_config(fasta, config_base, **kwargs)
#print ("Config commandline options:")
#print (config_stdout)
#
#cmd = (
#       "python {0}/ugly_hack_break_chain.py "
#       "--config {1}.up "
#       "--chain-break-from-file "
#      ).format(upside_utils_dir, pdb_id )
#print (cmd)
#sp.check_output(cmd.split())

