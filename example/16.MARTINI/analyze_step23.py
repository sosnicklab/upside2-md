import numpy as np
import tables as tb
import sys

def analyze_step23():
    """Analyze the trajectory to find which particles cause the +2.90 potential at step 23"""
    
    # MARTINI parameters
    epsilon = 1.715293655572
    sigma = 4.7
    
    # Open trajectory file
    traj_file = "./outputs/martini_test/test.run.up"
    
    with tb.open_file(traj_file, 'r') as f:
        # Read trajectory positions
        pos_data = f.root.output.pos[:]  # Shape: (n_frames, 1, n_atoms, 3)
        
        print(f"Trajectory shape: {pos_data.shape}")
        
        # Handle the extra dimension in UPSIDE output
        if len(pos_data.shape) == 4:
            pos_data = pos_data[:, 0, :, :]  # Remove the middle dimension
        
        n_frames, n_atoms, _ = pos_data.shape
        
        if n_frames <= 23:
            print(f"Not enough frames! Only {n_frames} frames available")
            return
            
        # Get positions at step 23 (frame 23)
        pos_23 = pos_data[23]  # Shape: (n_atoms, 3)
        
        print(f"\n=== ANALYSIS OF STEP 23 ===")
        print(f"Positions at step 23:")
        for i, pos in enumerate(pos_23):
            print(f"  Atom {i}: ({pos[0]:8.3f}, {pos[1]:8.3f}, {pos[2]:8.3f})")
        
        # Calculate all pairwise distances and LJ potentials
        total_potential = 0.0
        pair_contributions = []
        
        for i in range(n_atoms):
            for j in range(i+1, n_atoms):
                # Calculate distance
                dr = pos_23[i] - pos_23[j]
                dist = np.linalg.norm(dr)
                
                # Calculate LJ potential for this pair
                if dist > 0:
                    sig_r = sigma / dist
                    sig_r6 = sig_r**6
                    sig_r12 = sig_r6 * sig_r6
                    lj_pot = 4.0 * epsilon * (sig_r12 - sig_r6)
                else:
                    lj_pot = float('inf')
                
                total_potential += lj_pot
                pair_contributions.append((i, j, dist, lj_pot))
        
        print(f"\nTotal LJ potential: {total_potential:.3f}")
        
        # Sort by contribution (most positive first)
        pair_contributions.sort(key=lambda x: x[3], reverse=True)
        
        print(f"\n=== TOP 10 MOST REPULSIVE PAIRS ===")
        for i, (atom1, atom2, dist, pot) in enumerate(pair_contributions[:10]):
            print(f"{i+1:2d}. Atoms {atom1}-{atom2}: dist={dist:6.3f}Å, LJ_pot={pot:8.3f}")
        
        # Check for extremely close pairs
        print(f"\n=== CRITICALLY CLOSE PAIRS (dist < 2σ = {2*sigma:.1f}Å) ===")
        critical_pairs = [(a1, a2, d, p) for a1, a2, d, p in pair_contributions if d < 2*sigma]
        
        if critical_pairs:
            for atom1, atom2, dist, pot in critical_pairs:
                print(f"  Atoms {atom1}-{atom2}: dist={dist:.3f}Å, LJ_pot={pot:.3f}")
                print(f"    σ/r = {sigma/dist:.3f}, (σ/r)^6 = {(sigma/dist)**6:.3f}")
        else:
            print("  No critically close pairs found!")
        
        # Check minimum distance
        min_dist = min(d for _, _, d, _ in pair_contributions)
        max_pot = max(p for _, _, _, p in pair_contributions)
        
        print(f"\n=== SUMMARY ===")
        print(f"Minimum distance: {min_dist:.3f} Å")
        print(f"Maximum LJ potential: {max_pot:.3f}")
        print(f"σ (MARTINI): {sigma:.1f} Å")
        print(f"Minimum distance / σ: {min_dist/sigma:.3f}")

if __name__ == "__main__":
    analyze_step23() 