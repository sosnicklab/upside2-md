Restart an Upside Simulation Properly

# Procedures
1. first run `0.run.py` for 1000 steps
2. run `0.continue.py`, upside will correctly load the end configuration and the end momentum to properly restart the trajectory. (the only difference between `0.run.py` and `0.continue.py` is line 28 is changed to `continue_sim=True`).


# Explanation
1. in the first run, we call upside with the argument `--record-momentum` so that momenta recorded.
2. in the `0.continue.py`, line 28, with `continue_sim=True`, the python script will copy the configuration and the momentum of the end frame as the input position and momentum. 
3. in the `0.continue.py`, with `continue_sim=True`, the python script will call upside with the argument `--restart-using-momentum`

