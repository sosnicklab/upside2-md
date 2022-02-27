This example shows how to analyze Upside trajectories

1, use "mdtraj_upside" in python to import as a "Trajectory" of "mdtraj"
2, use the "tables" library in python to directly read the hdf5 file

After you finish the simulation, run 

    ./1.traj_ana.sh 

it will call the script "upside/py/get_info_from_upside_traj.py" to extract 
Rmsd, energy, Rg, #H-bond, Temperature. We encourage you to read this script 
to understand how we do that.

Then, you can compute any statistical property directly or with MBAR

    python 2.mbar_meltingCurve_freeEnergy.py


TODO:
add more comments in the scripts
add more explanation about MBAR
