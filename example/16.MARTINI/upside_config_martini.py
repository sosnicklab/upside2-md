# In your configuration script
kwargs = dict(
    overlayed_potential = {
        'epsilon': 1.0,
        'sigma': 1.0,
        'lj_cutoff': 10.0,
        'charges': [0.5, -0.5],  # Charges for each atom
        'coul_cutoff': 12.0,
        'dielectric': 80.0,
        'atom_indices': [0, 1]  # Indices of atoms to apply potential to
    }
)

# Use in simulation
config_stdout = ru.upside_config(fasta, config_base, 
    overlayed_potential=param_dir + "overlayed_potential.h5",
    **kwargs)
