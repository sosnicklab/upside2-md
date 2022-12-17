## Installation

### Dependencies

Compile dependencies

  * CMake, 2.8+
  * C++11 compiler, such as GCC 4.8+
  * HDF5, 1.8+ with high-level interface compile option
  * Eigen Matrix Library, 3.0+

Python dependencies.  Depending on your use-case for Upside, you may not need
all of these libraries.  The error message in case of a missing library should
be clear enough to figure out.

  * Python 3.7+
  * Numpy
  * Scipy (sometimes needed for `upside_config.py`, depending on options)
  * PyTables (needed for reading HDF5 files from Python)
  * Prody (needed for reading .pdb files)
  * Pandas (needed for `predict_chi1.py`)
  * H5py (needed for `generate_restart_config.py`)
  * MDTraj (optional but capable of loading and analyzing Upside trajectories)
  * Pymbar (optional but capable of doing mbar reweighting for REMD simulations)

### Compiling Upside

This guide assumes that you have all of the compile dependencies satisfied.  In
particular, the HDF5 high-level must be specifically enabled by a configuration
switch when installing HDF5.  Please check with your system administrator to
ensure it is installed.

This install guide assumes that you are using Linux or OS X.  Upside should
install and run on Windows as well, but this has not been tested.

The build system uses CMake (www.CMake.org) to find dependencies.  Please
change to the root directory of the upside directory. You may change the path 
for Eigen and python in source_sh based on your environment. Then, just run
the following commands.

    ./install.sh

After these commands execute successfully, the `obj/` directory will contain
the `upside` executable and the `libupside.so` shared library (exact name of
shared library may depend on operating system).

## Running molecular dynamics

The following sections illustrate running simple molecular dynamics simulations
with Upside.

### Converting a PDB file to Upside input

Starting from a PDB file, we create FASTA, position, and side chain files using

    upside/py/PDB_to_initial_structure.py /path/to/input.pdb output_basename

which will create `output_basename.fasta` containing the FASTA sequence,
`output_basename.initial.pkl` containing a Python pickle file of the
coordinates of the N, CA, and C atoms of the backbone, and
`output_basename.chi` containing the chi1 and chi2 angles for the side chains.
The .chi file is not needed for MD simulation, unless you want to run with
inflexible side chain rotamers.  `PDB_to_initial_structure.py` will fail if there
are breaks in the chain, as chain breaks will cause the MD to be invalid.  See
`--help` for information on how to choose the PDB model or chains to include.

From the output of `PDB_to_inital_structure.py`, we can create an Upside
configuration, `simulation.up`.  The following instructions are for a typical
folding or conformational dynamics simulation.

    upside/py/upside_config.py --output                   simulation.up \
                               --fasta                    output_basename.fasta \
                               --initial-structure        output_basename.initial.npy \
                               --hbond-energy             upside/parameters/ff_2.1/hbond.h5 \
                               --dynamic-rotamer-1body    \
                               --rotamer-placement        upside/parameters/ff_2.1/sidechain.h5 \
                               --rotamer-interaction      upside/parameters/ff_2.1/sidechain.h5 \
                               --environment-potential    upside/parameters/ff_2.1/environment.h5 \
                               --bb-environment-potential upside/parameters/ff_2.1/bb_env.dat \
                               --rama-library             upside/parameters/common/rama.dat \
                               --rama-sheet-mixing-energy upside/parameters/ff_2.1/sheet \
                               --reference-state-rama     upside/parameters/common/rama_reference.pkl

If `--initial-structure` is omitted, the simulation will be initialized with
random (possibly-clashing) Ramachandran angles, which is useful for de novo
structure prediction.

The output `simulation.up` is an HDF5 file, and the simulation configuration is
written in the "/input" group within the `.up` file.  The `upside_config.py`
copies all of the parameters into its output, so that only `simulation.up` is
needed to run the simulation.

### Constant temperature simulation

A simple, constant-temperature simulation may be run with 

    upside/obj/upside --duration 1e7 --frame-interval 1e2 --temperature 0.85 --seed $RANDOM simulation.up

Note that the simulation results are written into the "/output" group of the
provide `simulation.up` file.  This ensures that simulation output is stored in
the same file as the input, aiding in later reproducibility of the simulation.
The temperature, duration, and frame-interval are in natural units.  The times
should not be interpreted as picoseconds.  We are still investigating the
precise relationship of Upside time to experimental folding time, but the time
relationship may be conformation-dependent since the Upside backbone moves in a
smoother energy landscape than standard MD due to the side chain model. 

### Replica exchange simulation

The equilibration time of Upside simulations is highly temperature-dependent.
Typically, folding studies in Upside start with a replica exchange simulation
to establish the temperature range of the melting transition.  Hamiltonian and
temperature replica exchange are handled in the same manner in Upside.  First,
use `upside_config.py` repeatedly to create `N` configurations (or just copy
the output of a single `upside_config.py` run).  To run the replica exchange on
five configurations,

    upside/obj/upside --duration 1e7 --frame-interval 1e2 \
                      --temperature 0.5,0.53,0.56,0.60,0.64 \
                      --swap-set 0-1,2-3 \
                      --swap-set 1-2,3-4 \
                      --replica-interval 20 \
                      --monte-carlo-interval 5 \
                      --seed $RANDOM \
                      config_0.up config_1.up config_2.up config_3.up config_4.up

The multiple simulations are parallelized using OpenMP threads on a single
machine.  On the SLURM scheduler, jobs should be launched with `sbatch --ntasks
1 --cpus-per-task 5` for a five simulation replica exchange job.  The
`--replica-interval` controls the frequency at which replica exchange swaps are
attempted.  The `--swap-sets` arguments are *non-overlapping* pairs of replica
swaps that are attempted simultaneously.  For a linear chain of temperatures,
there should always be two swap sets to ensure ergodicity of swaps without
having any overlapping pairs.

The `--monte-carlo-interval` controls the frequency at which pivot moves, a
discrete change in the phi and psi angles of a single, randomly-chosen residue,
are attempted.  These Monte Carlo pivot moves are accepted with a Metropolis
criterion so that the Boltzmann ensemble of simulation is unchanged by pivot
moves.  Typically, pivot moves are highly advantageous in replica exchange
simulation because pivot moves have high acceptance rates in high-temperature
unfolded conformations, greatly speeding the sampling of unfolded states.
Pivot moves have very low acceptance in collapsed conformations but add very
little to the computational time.  The faster sampling in the unfolded states
feeds conformational diversity to the lower temperatures in replica exchange.
Monte Carlo moves are also permitted in constant temperature simulation.  If
`--monte-carlo-interval` is not provided, pivot moves are never attempted.

The output of replica exchange is written into the "/output" group of each of
the input configurations.  The replicas are written by temperature/Hamiltonian,
so that in the example above, `config_0.up` contains all of the simulation
frames at the temperature 0.5.  This means that in the limit of infinite sampling, 
`config_0.up` will contain samples from a Boltzmann ensemble at temperature 0.5 but 
the trajectory will be discontinuous due to replica swapping.

### Advanced usage

Explore the `upside/example/` folder for examples on running with restraints, membrane
simulations, pulling simulations, HDX prediction, and more.

## Simulation analysis and visualization

The contents of a `.up` file can be listed with

    upside/py/attr_overview.py simulation.up

which included copies of the full command lines used to invoke both
`upside_config.py` and `upside`.  

To load the simulation in VMD, first convert to a trajectory format
that VMD can read.  This can be done with

    upside/py/extract_vtf simulation.up simulation.vtf

and the `.vtf` file can be read natively by VMD. Alternatively, the MDTraj
library described below can read `simulation.up` directly, either to visualize
in an IPython notebook or to convert to another format.  In either case, amide
hydrogen and carbonyl oxygen are added as atoms in the trajectory to aid
structure viewers.

### Using MDTraj

To load the Upside trajectory as an MDTraj Trajectory object, use the
`mdtraj_upside.py` library in the Upside distribution.  

    import sys
    sys.path.append('upside/py')
    import mdtraj_upside as mu
    traj = mu.load_upside_traj('simulation.up')

When the trajectory is loaded, the atoms H, O, and CB will be added to the
trajectory.
