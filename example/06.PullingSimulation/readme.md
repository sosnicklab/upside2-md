This example shows how to perform pulling dynamics for membrane proteins.

Under this folder, just run the following commands to generate a trajectory:

    source ../../source.sh
    python run.py

then run

    python 1.get_force.py 1qhj_AFM.dat

to get the extension and force data to plot the extension force curve. The 
results are stored in 'results' folder.

There are two control files (1qhj_AFM.dat and 1qhj_AFM_auto_init.dat) for 
the pulling simulation. The difference between them is whether the initial 
tip_pos is specified. If not (such as 1qhj_AFM_auto_init.dat), the pulled 
atom position will be set to the initial tip_pos.


TODO:
