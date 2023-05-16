Restart an Upside Simulation Properly

# Procedures
1. first run `0.run.py` for 1000 steps
2. change the variable `continue_sim` from `False` to `True` in line 28.
3. run `0.run.py` again, upside will correctly load the end configuration and the end momentum to properly restart the trajectory. 

# Explanation
1. in the first run, we call upside with the argument `--record-momentum` so that momenta recorded.
2. in the `0.run.py`, line 204, if the `continue_sim=True`, the python script will copy the configuration and the momentum of the end frame as the input position and momentum. 
3. in the `0.run.py`, if the `continue_sim=True`, the python script will call upside with the argument `--restart-using-momentum`

