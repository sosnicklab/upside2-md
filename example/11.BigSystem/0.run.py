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

pdb_id         = 'S10000'
pdb_dir        = './pdb'
fasta_dir      = './fasta'

sim_id         = 'mem_test'
is_native      = False
ff             = 'ff_2.1'
T              = 0.8
duration       = 100
frame_interval = 10
base_dir       = './'

# Hi, HERE!, trj True, False respectively to see what happens
intensive_memory = False

n_rep            = 1     # replica number
randomseed       = 1     # np.random.randint(0,100000)
                         # Might want to change the fixed seed for the random number
account          = "your_account"    # FIXME change it 
partition        = "yout_partition"  # FIXME change it
#partition        = "caslake"
job_name         = '{}_{}'.format(pdb_id, sim_id)
run_time         = "36:00:00" # requested run time of job allocation in hh:mm:ss

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

# we turn off this step, just use .fasta file to make the config file
#print ("Initial structure gen...")
#cmd = (
#       "python {0}/PDB_to_initial_structure.py "
#       "{1}/{2}.pdb "
#       "{3}/{2} "
#       "--record-chain-breaks "
#      ).format(upside_utils_dir, pdb_dir, pdb_id, input_dir )
#print (cmd)
#sp.check_output(cmd.split())


#----------------------------------------------------------------------
## Configure
#----------------------------------------------------------------------

# parameters
param_dir_base = os.path.expanduser(upside_path+"/parameters/")
param_dir_common = param_dir_base + "common/"
param_dir_ff = param_dir_base + '{}/'.format(ff)

# options
print ("Configuring...")
fasta = "{}/{}.fasta".format(fasta_dir, pdb_id)
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
               intensive_memory          = intensive_memory,
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

h5_file  = "{}/{}.run.up".format(run_dir, pdb_id)
log_file = "{}/{}.run.log".format(run_dir, pdb_id)
shutil.copyfile(config_base, h5_file)

# SLURM options
# Will want to increase the time for production runs 
sbatch_opts = (
                "--account={} " #bphs35001
                "--job-name={} "
                "--output={} "
                "--time={} "
                "--partition={} "
                "--nodes=1 "
                "--ntasks-per-node={} "
                "--mem=8000"
              )

sbatch_opts = sbatch_opts.format(account, job_name, log_file, run_time, partition, n_rep)

print ("Running...")
cmd = "sbatch {} --wrap=\"{}/obj/upside {} {}\"".format(sbatch_opts, upside_path, upside_opts, h5_file)
sp.check_call(cmd, shell=True)
