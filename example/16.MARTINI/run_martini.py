import sys, os, shutil
import subprocess as sp
import numpy as np
import tables as tb

# Get UPSIDE home directory
upside_path = os.environ['UPSIDE_HOME']
upside_utils_dir = os.path.expanduser(upside_path+"/py")
sys.path.insert(0, upside_utils_dir)

# Get the Python interpreter path
python_path = sys.executable

# Import UPSIDE utilities
import run_upside as ru
from advanced_config import write_external_pairs_potential

#----------------------------------------------------------------------
## General Settings and Path
#----------------------------------------------------------------------

sim_id         = 'martini_test'
T              = 0.8  # Reduced temperature
duration       = 1000
frame_interval = 50
base_dir       = './'

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
## Configure
#----------------------------------------------------------------------

# parameters
param_dir_base = os.path.expanduser(upside_path+"/parameters/")
param_dir_common = param_dir_base + "common/"
param_dir_ff = param_dir_base + 'ff_2.1/'  # Using standard force field

# options
print("Configuring...")
fasta = "{}/martini.fasta".format(input_dir)  # You'll need to create this

# Create a minimal fasta file with the correct number of residues
with open(fasta, 'w') as f:
    f.write("> MARTINI\n")
    # Add appropriate number of 'A' residues based on your system
    f.write("A" * 100)  # Adjust number as needed

# Build command line arguments for upside_config
config_args = [
    python_path,
    os.path.join(upside_utils_dir, 'upside_config.py'),
    '--fasta', fasta,
    '--output', "{}/martini.up".format(input_dir),
    '--rama-library', param_dir_common + "rama.dat",
    '--rama-sheet-mixing-energy', param_dir_ff + "sheet",
    '--reference-state-rama', param_dir_common + "rama_reference.pkl",
    '--hbond-energy', param_dir_ff + "hbond.h5",
    '--rotamer-placement', param_dir_ff + "sidechain.h5",
    '--dynamic-rotamer-1body',
    '--rotamer-interaction', param_dir_ff + "sidechain.h5",
    '--environment-potential', param_dir_ff + "environment.h5",
    '--bb-environment-potential', param_dir_ff + "bb_env.dat",
]

print("Config commandline options:")
print(' '.join(config_args[2:]))  # Skip python path and script name

# Run upside_config
try:
    sp.check_call(config_args)
except sp.CalledProcessError as e:
    print(f"Error running upside_config: {e}")
    sys.exit(1)

# Add custom potential using write_external_pairs_potential
with tb.open_file("{}/martini.up".format(input_dir), 'a') as t:
    # Set the global variable t that write_external_pairs_potential expects
    import advanced_config
    advanced_config.t = t
    
    # Set up the potential group
    if not hasattr(t.root.input, 'potential'):
        t.create_group(t.root.input, 'potential')
    advanced_config.potential = t.root.input.potential
    
    # First write CB positions
    with open(fasta, 'r') as f:
        fasta_content = f.read().split('\n')[1]  # Skip header line
    advanced_config.write_CB(fasta_content)
    
    # Create overlayed_potential group
    grp = t.create_group(t.root.input.potential, 'overlayed_potential')
    grp._v_attrs.arguments = np.array([b'pos'])
    
    # Read data from MARTINI h5 file
    with tb.open_file("{}/martini.h5".format(input_dir)) as epp:
        pairs = epp.root.pairs[:]
        coefficients = epp.root.coefficients[:]
        atoms = epp.root.atoms[:]
        
        # Convert atom types
        atoms[atoms=='CB'] = 0
        atoms[atoms=='CA'] = 1
        atoms = np.array(atoms, dtype=int)
        atoms = atoms.flatten()  # Flatten to 1D array
        
        # Set LJ and Coulombic parameters
        grp._v_attrs.epsilon = 1.0  # LJ epsilon
        grp._v_attrs.sigma = 1.0    # LJ sigma
        grp._v_attrs.lj_cutoff = 1.2  # LJ cutoff distance
        grp._v_attrs.coul_cutoff = 1.2  # Coulombic cutoff distance
        grp._v_attrs.dielectric = 1.0  # Dielectric constant
        
        # Create arrays in the overlayed_potential group
        tb.Array(grp, 'atom_indices', obj=atoms)
        tb.Array(grp, 'charges', obj=np.zeros(len(atoms)))  # Zero charges for now

#----------------------------------------------------------------------
## Run Settings
#----------------------------------------------------------------------

randomseed = 1  # Fixed seed for reproducibility

upside_opts = (
    "--duration {} "
    "--frame-interval {} "
    "--temperature {} "
    "--seed {} "
)
upside_opts = upside_opts.format(duration, frame_interval, T, randomseed)

h5_file  = "{}/martini.run.up".format(run_dir)
log_file = "{}/martini.run.log".format(run_dir)
shutil.copyfile("{}/martini.up".format(input_dir), h5_file)

print("Running...")

# Check if we're on a SLURM system
is_slurm = 'SLURM_JOB_ID' in os.environ

if is_slurm:
    # SLURM system - use srun
    cmd = "srun {}/obj/upside {} {} | tee {}".format(upside_path, upside_opts, h5_file, log_file)
else:
    # Non-SLURM system (like Mac) - run directly
    cmd = "{}/obj/upside {} {} | tee {}".format(upside_path, upside_opts, h5_file, log_file)

# Execute the command
try:
    sp.check_call(cmd, shell=True)
except sp.CalledProcessError as e:
    print(f"Error running UPSIDE: {e}")
    sys.exit(1)
except FileNotFoundError as e:
    print(f"Error: Could not find UPSIDE executable. Make sure UPSIDE_HOME is set correctly.")
    print(f"Current UPSIDE_HOME: {upside_path}")
    sys.exit(1)
