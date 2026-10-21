Self Avoid Random Walk Simulation

In this simulation, only the backbone terms are kept: the bond, 
angle and omega springs, the Ramachandran potential, and the 
steric repulsion between backbone atoms, which makes the walk 
self-avoiding. The non-bonded terms (hydrogen bonds, side chains 
and burial environment) are left out. Since all the expensive 
forces are excluded, the simulation will be much faster than a 
normal simulation.

TODO:
add the script to calculate the P_of_r
