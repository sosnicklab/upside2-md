This example shows how to calculate HDX and its denaturant dependencies 
based on REMD trajectories of Upside
It is best to try example 3 first.

step 1. run REMD simulations 

    source ../../source.sh
    python 0.run.py

 This step is exactly the same as the first step in Example 3, 
 so you can copy the inputs and outputs of 3.TrajectoryAnalysis 
 to here directly.


step 2. configure a new .up file (inputs/EHEE_rd2_0005-HDX.up) which
        allows to calculate the burial level of backbone NH

    python 1.config.py

step 3. get T, energy in Upside trajs to build MBER

    sh 2.traj_ana.sh

step 4. get the protaction states for all residues in evergy frame

    sh 3.get_protaction_states.sh

step 5. calculate the dG_HX and the m-value

    python 4.calc_HDX.py


TODO:
