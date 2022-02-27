import sys, os, time, shutil
from math import sqrt
import subprocess as sp
import numpy as np
import tables as tb


upside_path = '/home/pengxd/upside2/'
upside_utils_dir = os.path.expanduser(upside_path+"/py")
sys.path.insert(0, upside_utils_dir)
import run_upside as ru
from scale_params import *

## General Settings and Path
pdb_id    = sys.argv[1]
su        = sys.argv[2]
is_native = True

# Parameters
param_dir_base = os.path.expanduser(upside_path+"/parameters/")
param_dir_common = param_dir_base + "common/"
param_dir_ff1 = param_dir_base + 'ff_2.0/'

## Generate Upside readable initial structure (and fasta) from PDB 
print ("Initial structure gen...")
cmd = (
       "python {0}/PDB_to_initial_structure.py "
       "pdb/{1}.pdb "
       "{1} "
       "--record-chain-breaks "
       "--allow-unexpected-chain-breaks "
       "--disable-recentering "
      ).format(upside_utils_dir, pdb_id )
print (cmd)
sp.check_output(cmd.split())

## Configure
print ("Configuring...")
fasta = "{}.fasta".format(pdb_id)
kwargs = dict(
               rama_pot                  = param_dir_common + "rama.dat",
               hbond                     = param_dir_ff1 + "hbond.h5",
               sheet_mix_energy          = param_dir_ff1 + "sheet",
               reference_rama            = param_dir_common + "rama_reference.pkl",
               placement                 = param_dir_ff1 + "sidechain.h5",
               dynamic_1body             = True,
               rotamer_interaction_param = param_dir_ff1 + "sidechain.h5",
               environment               = param_dir_ff1 + "environment.h5",
               environment_type          = 1,
               chain_break_from_file     = "{}.chain_breaks".format(pdb_id),
               bb_environment            = param_dir_ff1 + "bb_env.dat",
               #trans_cis                 = 'trans_cis.dat',
             )

if is_native:
    kwargs['init'] =  "{}.initial.npy".format(pdb_id)

config_base = "{}-{}.up".format( pdb_id, su)
config_stdout = ru.upside_config(fasta, config_base, **kwargs)

print ("Config commandline options:")
print (config_stdout)

## advanced configure
config_base = "{}-{}.up".format( pdb_id, su)
kwargs = dict(
               #restraint_groups = ['0-30','31-60'],
               AFM = 'pulling.dat',
               #fixed_wall = 'wall-{}.dat'.format(su),
               #pair_wall  = 'wall-{}.dat'.format(su),
               #fixed_spring  = 'spring-{}.dat'.format(su),
               #nail  = 'nail.dat',
             )

config_stdout = ru.advanced_config(config_base, **kwargs)
print ("Advanced Config commandline options:")
print (config_stdout)
