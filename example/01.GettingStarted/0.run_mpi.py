import sys, os, shutil
import subprocess as sp
import numpy as np

upside_path = os.environ['UPSIDE_HOME']
upside_utils_dir = os.path.expanduser(upside_path+"/py")
sys.path.insert(0, upside_utils_dir)
import run_upside as ru

#----------------------------------------------------------------------
## General Settings and Path
#----------------------------------------------------------------------

pdb_id         = 'chig' # switch to 1dfn for multi-chain showing
pdb_dir        = './pdb'
sim_id         = 'simple_test_mpi'
is_native      = True
ff             = 'ff_2.1'
T              = 0.8
duration       = 1000
frame_interval = 50
base_dir       = './'

# MPI settings
use_mpi = True  # Set to True to use MPI
mpi_ranks = 2   # Number of MPI ranks to use

#----------------------------------------------------------------------
## Initialization
#----------------------------------------------------------------------

input_dir  = "{}/inputs".format(base_dir)
output_dir = "{}/outputs".format(base_dir)
run_dir    = "{}/{}".format(output_dir, sim_id) 

make_dirs = [input_dir, output_dir, run_dir]
for direc in make_dirs:
    if not os.path.exists(direc):
        os.makedirs(direc)

#----------------------------------------------------------------------
## Generate Upside readable initial structure (and fasta) from PDB 
#----------------------------------------------------------------------

print ("Initial structure gen...")
cmd = (
       "python {0}/PDB_to_initial_structure.py "
       "{1}/{2}.pdb "
       "{3}/{2} "
       "--record-chain-breaks "
      ).format(upside_utils_dir, pdb_dir, pdb_id, input_dir )
print (cmd)
sp.check_output(cmd.split())


#----------------------------------------------------------------------
## Configure
#----------------------------------------------------------------------

# parameters
param_dir_base = os.path.expanduser(upside_path+"/parameters/")
param_dir_common = param_dir_base + "common/"
param_dir_ff = param_dir_base + '{}/'.format(ff)

# options
print ("Configuring...")
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
             )

if is_native:
    kwargs['initial_structure'] =  "{}/{}.initial.npy".format(input_dir, pdb_id)

config_base = "{}/{}.up".format( input_dir, pdb_id)
config_stdout = ru.upside_config(fasta, config_base, **kwargs)

print ("Config commandline options:")
print (config_stdout)

#----------------------------------------------------------------------
## Run Settings
#----------------------------------------------------------------------

randomseed = np.random.randint(0,100000) # Might want to change the fixed seed for the random number
randomseed = 1

upside_opts = (
                 "--duration {} "
                 "--frame-interval {} "
                 "--temperature {} "
                 "--seed {} "
              )
upside_opts = upside_opts.format(duration, frame_interval, T, randomseed)

# MPI setup
if use_mpi:
    print(f"Using MPI with {mpi_ranks} ranks")
    mpi_cmd = f"mpirun -np {mpi_ranks}"
    
    # Create multiple copies of the system for MPI parallelization
    system_files = []
    for i in range(mpi_ranks):
        copy_file = "{}/{}.run.{}.up".format(run_dir, pdb_id, i)
        shutil.copyfile(config_base, copy_file)
        system_files.append(copy_file)
    
    h5_file = " ".join(system_files)
    print(f"Created {mpi_ranks} system copies for MPI parallelization")
else:
    print("Using single-process execution")
    mpi_cmd = ""
    h5_file = "{}/{}.run.up".format(run_dir, pdb_id)
    shutil.copyfile(config_base, h5_file)

log_file = "{}/{}.run.log".format(run_dir, pdb_id)

print ("Running...")
if use_mpi:
    cmd = "{} {}/obj/upside {} {} | tee {}".format(mpi_cmd, upside_path, upside_opts, h5_file, log_file)
else:
    cmd = "{}/obj/upside {} {} | tee {}".format(upside_path, upside_opts, h5_file, log_file)

print(f"Command: {cmd}")
sp.check_call(cmd, shell=True)

print("MPI simulation completed successfully!")
