import sys, os, shutil
import subprocess as sp
import numpy as np
import tables as tb
from math import sqrt
import time

upside_path = os.environ['UPSIDE_HOME']
upside_utils_dir = os.path.expanduser(upside_path+"/py")
sys.path.insert(0, upside_utils_dir)
import run_upside as ru

#----------------------------------------------------------------------
## General Settings and Path
#----------------------------------------------------------------------

pdb_id           = 'three_1rkl_a'
pdb_dir          = './pdb'
sim_id           = 'curv_dy1'
is_native        = True
ff               = 'ff_2.1'
thickness        = 24.8

duration         = 100001
frame_interval   = 100
work_dir         = './'

curvature_radius = 250.  # curvature radius (angstrom)
curvature_sign   = 1     # 1 means positive curvature, -1 means negative curvature

exchange         = False # if True, it will run the replica exchange simulation
                         # if False, it will run the constant temperature simulation

n_rep            = 4     # replica number
T_low            = 0.80 
T_high           = 0.80
replica_interval = 10    # How long takes an exchange attempt (upside time unit)

continue_sim     = False # when you run a new simulation, set it as "False"
                         # "True" means restarting the simulation from the last frame
                         # of the previous trajectories (they should have the same 
                         # pdb_id and sim_id as the new simulation, and exist in the 
                         # corresponding path)

randomseed       = -1    # np.random.randint(0,100000) 
                         # Might want to change the fixed seed for the random number

account          = "your_account"    # FIXME change it 
partition        = "yout_partition"  # FIXME change it
job_name         = '{}_{}'.format(pdb_id, sim_id)
run_time         = "36:00:00" # requested run time of job allocation in hh:mm:ss


#----------------------------------------------------------------------
# Set the path and filename
#----------------------------------------------------------------------

input_dir  = "{}/inputs".format(work_dir)
output_dir = "{}/outputs".format(work_dir)
run_dir    = "{}/{}".format(output_dir, sim_id) 

make_dirs = [input_dir, output_dir, run_dir]
for direc in make_dirs:
    if not os.path.exists(direc):
        os.makedirs(direc)

h5_files = []
for j in range(n_rep): 
    h5_file  = "{}/{}.run.{}.up".format(run_dir, pdb_id, j)
    h5_files.append(h5_file)
h5_files_str = " ".join(h5 for h5 in h5_files)
log_file = "{}/{}.run.log".format(run_dir, pdb_id)

#----------------------------------------------------------------------
## Check the previous trajectories if you set continue_sim = True 
#----------------------------------------------------------------------

if continue_sim:
    for h5 in h5_files:
        exist = os.path.exists(h5)
        if not exist:
            print('Warning: no previous trajectory file {}!'.format(h5))
            print('set "continue_sim = False" and start a new simulation')
            continue_sim = False
            break
    if continue_sim:
        exist = os.path.exists(log_file)
        if not exist:
            print('Warning: no previous log file {}!'.format(log_file))

#----------------------------------------------------------------------
## Generate Upside readable initial structure (and fasta) from PDB 
#----------------------------------------------------------------------

if not continue_sim:
    print ("Initial structure gen...")
    cmd = (
           "python {0}/PDB_to_initial_structure.py "
           "{1}/{2}.pdb "
           "{3}/{2} "
           "--record-chain-breaks "
           "--disable-recentering "
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
               membrane_potential        = param_dir_ff + "membrane.h5",
               membrane_thickness        = thickness,
               chain_break_from_file     = "{}/{}.chain_breaks".format(input_dir, pdb_id),
               use_curvature             = True,
               curvature_radius          = curvature_radius,
               curvature_sign            = curvature_sign,
             )

if is_native:
    kwargs['initial_structure'] =  "{}/{}.initial.npy".format(input_dir, pdb_id)

config_base = "{}/{}.up".format( input_dir, pdb_id)

if not continue_sim:
    print ("Configuring...")
    config_stdout = ru.upside_config(fasta, config_base, **kwargs)
    print ("Config commandline options:")
    print (config_stdout)

if not continue_sim:
    kwargs = dict(
                   nail = 'nail-xyz.dat'
                 )
    
    config_stdout = ru.advanced_config(config_base, **kwargs)
    print ("Advanced Config commandline options:")
    print (config_stdout)


#----------------------------------------------------------------------
## Run Settings
#----------------------------------------------------------------------

upside_opts = (
                 "--duration {} "
                 "--frame-interval {} "
                 "--temperature {} "
                 "--seed {} "
                 "--disable-recentering "
                 "--integrator mv "
                 "--curvature-changer-interval 10 "
               )

tempers =  np.linspace(sqrt(T_low), sqrt(T_high), n_rep)**2
tempers_str = ",".join(str(t) for t in tempers)

if exchange:
    swap_sets    = ru.swap_table2d(1, len(tempers)) # specifies which replicas are able to exchange 
    upside_opts += "--replica-interval {} --swap-set {} --swap-set {} " # only perform swaps for replex; duration of time until swap is attempted
    upside_opts  = upside_opts.format(duration, frame_interval, tempers_str, randomseed, replica_interval, swap_sets[0], swap_sets[1])
else:
    upside_opts  = upside_opts.format(duration, frame_interval, tempers_str, randomseed)


if continue_sim:
    print ("Archiving prev output...")

    localtime = time.asctime( time.localtime(time.time()) )
    localtime = localtime.replace('  ', ' ')
    localtime = localtime.replace(' ', '_')
    localtime = localtime.replace(':', '-')

    if os.path.exists(log_file):
        shutil.move(log_file, '{}.bck_{}'.format(log_file, localtime))
    else:
        print('Warning: no previous log file {}!'.format(log_file))

    for fn in h5_files:
        with tb.open_file(fn, 'a') as t:
            i = 0
            while 'output_previous_%i'%i in t.root:
                i += 1
            new_name = 'output_previous_%i'%i
            if 'output' in t.root:
                n = t.root.output
            else:
                n = t.get_node('/output_previous_%i'%(i-1))

            t.root.input.pos[:,:,0] = n.pos[-1,0]

            if 'output' in t.root:
                t.root.output._f_rename(new_name)

else:
    for fn in h5_files:
        shutil.copyfile(config_base, fn)

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
              )

sbatch_opts = sbatch_opts.format(account, job_name, log_file, run_time, partition, n_rep)

print ("Running...")
#cmd = "{}/obj/upside {} {}".format( upside_path, upside_opts, h5_files_str)
cmd = "sbatch {} --wrap=\"{}/obj/upside {} {}\"".format(sbatch_opts, upside_path, upside_opts, h5_files_str)
sp.check_call(cmd, shell=True)
