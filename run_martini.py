# Create martini_potential group
martini_group = t.create_group(t.root.input.potential, 'pair_potential')  # Use standard pair_potential name
martini_group._v_attrs.arguments = np.array([b'pos'])
martini_group._v_attrs.potential_type = b'pair'  # Use standard pair potential type
martini_group._v_attrs.epsilon = 1.715293655572  # Original MARTINI epsilon (kJ/mol)
martini_group._v_attrs.sigma = 4.7    # Original MARTINI sigma (Å)
martini_group._v_attrs.lj_cutoff = 12.0  # Å (converted from 1.2 nm)
martini_group._v_attrs.coul_cutoff = 12.0  # Å (converted from 1.2 nm)
martini_group._v_attrs.dielectric = 15.0  # MARTINI water dielectric
martini_group._v_attrs.n_types = 1  # Only one type of particle (water)
martini_group._v_attrs.n_params = 4  # [epsilon, sigma, charge1, charge2]
martini_group._v_attrs.cutoff = 12.0  # Å (converted from 1.2 nm)
martini_group._v_attrs.cache_buffer = 1.0  # Default cache buffer
martini_group._v_attrs.initialized = True  # Mark as initialized

# Remove all other potential groups except pair_potential
groups_to_remove = [group for group in t.root.input.potential._v_groups 
                   if group not in ['pair_potential']]
for group in groups_to_remove:
    t.remove_node(t.root.input.potential, group, recursive=True)

# Now copy the temporary file to the final location
import shutil
shutil.copy2(temp_file, "{}/martini.up".format(input_dir))

# Open the final file for modification
with tb.open_file("{}/martini.up".format(input_dir), 'r+') as t:
    # Create martini_potential group
    martini_group = t.create_group(t.root.input.potential, 'martini_potential')
    martini_group._v_attrs.arguments = np.array([b'pos'])
    martini_group._v_attrs.potential_type = b'pair'  # Use standard pair potential type
    martini_group._v_attrs.epsilon = 1.715293655572  # Original MARTINI epsilon (kJ/mol)
    martini_group._v_attrs.sigma = 4.7    # Original MARTINI sigma (Å)
    martini_group._v_attrs.lj_cutoff = 12.0  # Å (converted from 1.2 nm)
    martini_group._v_attrs.coul_cutoff = 12.0  # Å (converted from 1.2 nm)
    martini_group._v_attrs.dielectric = 15.0  # MARTINI water dielectric
    martini_group._v_attrs.n_types = 1  # Only one type of particle (water)
    martini_group._v_attrs.n_params = 4  # [epsilon, sigma, charge1, charge2]
    martini_group._v_attrs.cutoff = 12.0  # Å (converted from 1.2 nm)
    martini_group._v_attrs.cache_buffer = 1.0  # Default cache buffer
    martini_group._v_attrs.initialized = True  # Mark as initialized
    
    # Copy arrays from pair_potential to martini_potential
    if hasattr(t.root.input.potential, 'pair_potential'):
        pair_group = t.root.input.potential.pair_potential
        # Copy atom_indices
        t.create_array(martini_group, 'atom_indices', obj=pair_group.atom_indices[:])
        # Copy charges
        t.create_array(martini_group, 'charges', obj=pair_group.charges[:])
        # Copy pairs
        t.create_array(martini_group, 'pairs', obj=pair_group.pairs[:])
        # Copy coefficients
        t.create_array(martini_group, 'coefficients', obj=pair_group.coefficients[:])
        
        # Remove pair_potential group
        t.remove_node(t.root.input.potential, 'pair_potential', recursive=True) 