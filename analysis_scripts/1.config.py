import sys, os, shutil
import runpy
import subprocess as sp
import numpy as np
import tables as tb
from math import sqrt
import time

workflow_action = os.environ.get('workflow_action', 'config').strip().lower()
if workflow_action == 'select_h5_frames':
    runpy.run_module('helpers.select_h5_frames', run_name='__main__')
    raise SystemExit(0)
if workflow_action != 'config':
    raise SystemExit("Unsupported workflow_action for 1.config.py: {}".format(workflow_action))

upside_path = os.environ['UPSIDE_HOME']
upside_utils_dir = os.path.expanduser(upside_path+"/py")
sys.path.insert(0, upside_utils_dir)
import run_upside as ru

#----------------------------------------------------------------------
## General Settings and Path
#----------------------------------------------------------------------

pdb_id           = os.environ.get('pdb_id', 'glpG-RKRK-79HIS')  # CHECKME
pdb_dir          = './pdb'
is_native_env    = os.environ.get('is_native')
is_native        = True if is_native_env is None else is_native_env.lower() in ('1', 'true', 'yes', 'y')  # CHECKME
ff               = os.environ.get('ff', 'ff_2.1').strip()  # CHECKME
work_dir         = './'
input_dir        = "{}/inputs".format(work_dir)

#----------------------------------------------------------------------
## Configure
#----------------------------------------------------------------------

# parameters
param_dir_base = os.path.expanduser(upside_path+"/parameters/")
param_dir_common = param_dir_base + "common/"
param_dir_ff = param_dir_base + '{}/'.format(ff)

# options
fasta = "{}/{}.fasta".format(input_dir, pdb_id)
kwargs = dict(
               rama_library              = param_dir_common + "rama.dat",
               rama_sheet_mix_energy     = param_dir_ff + "sheet",
               reference_state_rama      = param_dir_common + "rama_reference.pkl",
               hbond_energy              = param_dir_ff + "hbond.h5",
               rotamer_placement         = param_dir_ff + "sidechain.h5",
               dynamic_rotamer_1body     = True,
               rotamer_interaction       = param_dir_ff + "sidechain.h5",
               environment_potential     = param_dir_ff + "environment.h5",
               bb_environment_potential  = param_dir_ff + "bb_env.dat",
               chain_break_from_file     = "{}/{}.chain_breaks".format(input_dir, pdb_id),
               use_heavy_atom_coverage   = True
             )

if is_native:
    kwargs['initial_structure'] =  "{}/{}.initial.npy".format(input_dir, pdb_id)

config_base = "{}/{}-HDX.up".format( input_dir, pdb_id)

print ("Configuring...")
config_stdout = ru.upside_config(fasta, config_base, **kwargs)
print ("Config commandline options:")
print (config_stdout)
