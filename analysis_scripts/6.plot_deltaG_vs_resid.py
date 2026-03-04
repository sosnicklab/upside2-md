import numpy as np
import matplotlib.pyplot as plt
import pymbar # for MBAR analysis
from pymbar import timeseries # for timeseries analysis
import subprocess
import os

pdb_id      = 'glpG-RKRK-79HIS'  #CHECKME
sim_id      = 'memb_test'  #CHECKME
start_frame = 100

work_dir         = './'
n_rep            = 48     # replica number

input_dir  = "{}/inputs".format(work_dir)
output_dir = "{}/outputs".format(work_dir)
result_dir = "{}/results".format(work_dir)
run_dir    = "{}/{}".format(output_dir, sim_id)

#============================================
# load data
#============================================

Pot = []
Rg  = []
Rmsd = []
Hb = []
Ts = []
T  = []
PS = []
for i in range(n_rep):
    Pot.append(  np.load('{}/{}_{}_{}_Energy.npy'.format(result_dir, pdb_id, sim_id, i))[:,0] ) 
    Rg.append(   np.load('{}/{}_{}_{}_Rg.npy'    .format(result_dir, pdb_id, sim_id, i)) ) 
    Hb.append(   np.load('{}/{}_{}_{}_Hbond.npy' .format(result_dir, pdb_id, sim_id, i)) ) 
    Rmsd.append( np.load('{}/{}_{}_{}_Rmsd.npy'  .format(result_dir, pdb_id, sim_id, i)) ) 
    PS.append( np.load('{}/{}_{}_{}_PS.npy'  .format(result_dir, pdb_id, sim_id, i)) ) 

    t = np.load( '{}/{}_{}_{}_T.npy'.format(result_dir, pdb_id, sim_id, i) )
    nsize = Pot[-1].size
    Ts.append(np.zeros(nsize) + t )
    T.append(t)
    
Rmsd = np.array(Rmsd)
Pot  = np.array(Pot)
Rg   = np.array(Rg)
Hb   = np.array(Hb)
Ts   = np.array(Ts)
T    = np.array(T)
PS   = np.array(PS)

res = np.loadtxt('{}/{}.resid'.format(result_dir, pdb_id), dtype=int)
n_res = res.size

print(PS.shape)

#============================================
# mbar
#============================================

kB   = 1.0 # upside unit
T    = np.array(T)
beta = kB*T**(-1)

cE0 = Pot[:,start_frame:]

FN           = cE0[0].size
FNs          = np.zeros([n_rep], np.int32) + FN
reducedPot0  = np.zeros([n_rep,n_rep,FN], np.float32)
for k in range(n_rep):
    for l in range(n_rep):
        reducedPot0[k,l] = beta[l] * cE0[k]
mbar0 = pymbar.MBAR(reducedPot0, FNs, verbose=True)

#============================================
# calculate dG_HX at different T 
#============================================
# (Keeping original logic for the basic dG array, though not used in final plot)
dGhx_T = []
for k in range(T.size):
    t = T[k]
    tt = t/0.85*298.

    u_n = (cE0/(t*kB)).flatten()
    log_w1 = mbar0._computeUnnormalizedLogWeights(u_n)
    w1 = np.exp(log_w1)
    w1 /= np.sum(w1)

    dG = np.zeros(n_res)
    for r in range(n_res):
        pf_i = PS[:,start_frame:,r].flatten()
        mean_pf = np.average(pf_i, weights=w1)
        if mean_pf == 1:
            # print(k, r)
            dG[r] = 1000.
        else:
            dG[r] = 0.001987*tt*np.log((mean_pf/(1.-mean_pf)))
    dGhx_T.append(dG)
dGhx_T = np.array(dGhx_T)

#============================================
# calculate dG_HX at different [den]
#============================================

T_target = 0.8   # the T you want to use 
m_sens   = 0.04  # it is the sensitivity of the denaturant.
                 # a larger value means a lower concentration 
                 # to unfold the protein
                 # 0.05 is a not bad initial guess for urea

# the number of protected residues in evergy frame
pf_frame = np.sum(PS[:,start_frame:,:], axis=2).flatten()

# [den]
den_min = 0.0
den_max = 15
den_bin = 150
den  = np.linspace(den_min, den_max, den_bin+1)

# reweight to T_target
TT = T_target/0.85*298
m = m_sens*-pf_frame
u_n    = (cE0/(T_target*kB)).flatten()
log_w1 = mbar0._computeUnnormalizedLogWeights(u_n)
w1     = np.exp(log_w1)
w1    /= np.sum(w1)

# dG at different [den]
dGhx_D = []
for d in den:
    w  = np.exp(m*d/T_target)*w1
    w /= np.sum(w)
    
    probp = np.zeros(n_res)
    for j in range(n_res):
        pf_i = PS[:,start_frame:,j].flatten()
        probp[j] = np.sum(pf_i*w)
    probu = 1. - probp
    dGhx_D.append(TT*np.log(probp/probu)*0.001987)
dGhx_D = np.array(dGhx_D)

# m-value: the slope of dG vs [den]
mValue = np.diff(dGhx_D, axis=0)/(den[0]-den[1])

#============================================
# plot residue ID vs delta G for different temperatures
# WITH ERROR ANALYSIS (Jackknife over Trajectories)
#============================================

# Target temperatures
target_temps = [0.70, 0.75, 0.8, 0.85, 0.90]

# Lists to store final mean and error for plotting
dGhx_target_mean = []
dGhx_target_err  = []

print("\nStarting Error Analysis (Jackknife over {} trajectories)...".format(n_rep))

for target_temp in target_temps:
    tt = target_temp/0.85*298.
    
    # 1. Get global MBAR weights for this temperature
    u_n = (cE0/(target_temp*kB)).flatten()
    log_w1 = mbar0._computeUnnormalizedLogWeights(u_n)
    w1 = np.exp(log_w1)
    
    # NOTE: w1 is shape (n_rep * n_frames). 
    # To compute error between trajectories, we need to separate contributions by trajectory.
    # The PS array is (n_rep, n_frames, n_res).
    
    # Reshape weights to match (n_rep, n_frames)
    frames_per_rep = cE0.shape[1]
    w1_matrix = w1.reshape(n_rep, frames_per_rep)
    
    # Pre-calculate sums for the ratio estimator
    # Denom for weighted avg: Sum of weights per trajectory
    sum_w_per_traj = np.sum(w1_matrix, axis=1) 
    
    dG_mean = np.zeros(n_res)
    dG_err  = np.zeros(n_res)
    
    for r in range(n_res):
        # Protection status for residue r: (n_rep, n_frames)
        pf_r = PS[:, start_frame:, r]
        
        # Numerator for weighted avg: Sum of (w * p) per trajectory
        sum_wp_per_traj = np.sum(pf_r * w1_matrix, axis=1)
        
        # --- Standard Calculation (All Data) ---
        total_num = np.sum(sum_wp_per_traj)
        total_den = np.sum(sum_w_per_traj)
        
        # Calculate protection fraction
        pf_all = total_num / total_den
        
        # Calculate dG (Mean)
        if pf_all >= 0.99999:
            dG_mean[r] = 1000.
            dG_err[r]  = 0.
            continue
        elif pf_all <= 0.00001:
            dG_mean[r] = -100. # Arbitrary low value for fully unfolded
            dG_err[r]  = 0.
            continue
        else:
            dG_mean[r] = 0.001987 * tt * np.log(pf_all / (1. - pf_all))

        # --- Jackknife Error Estimation ---
        # We recalculate dG n_rep times, removing one trajectory each time.
        jk_dG_values = []
        
        for k in range(n_rep):
            # Subtract the contribution of trajectory k
            jk_num = total_num - sum_wp_per_traj[k]
            jk_den = total_den - sum_w_per_traj[k]
            
            jk_pf = jk_num / jk_den
            
            # Guard against log errors in jackknife samples
            jk_pf = max(min(jk_pf, 1.0 - 1e-6), 1e-6)
            
            jk_val = 0.001987 * tt * np.log(jk_pf / (1. - jk_pf))
            jk_dG_values.append(jk_val)
            
        jk_dG_values = np.array(jk_dG_values)
        
        # Jackknife estimate of standard error
        # Formula: sqrt( (N-1)/N * sum( (theta_i - theta_mean)^2 ) )
        # where theta_i are the leave-one-out estimates.
        # np.var calculates mean((x-mean)^2), so we multiply by (N-1) to get the jackknife sum
        variance = (n_rep - 1) * np.var(jk_dG_values)
        dG_err[r] = np.sqrt(variance)

    dGhx_target_mean.append(dG_mean)
    dGhx_target_err.append(dG_err)
    print(f"Calculated error bars for T={target_temp}")

# Plot Delta G vs Residue ID with Error Bars
fig = plt.figure(figsize=(12, 8))
cmap = plt.cm.get_cmap('tab10', len(target_temps))

for i, (temp, dG_values, dG_errors) in enumerate(zip(target_temps, dGhx_target_mean, dGhx_target_err)):
    
    # Filter out the placeholder '1000' values for plotting clarity if desired, 
    # or plot them (they will be off scale). 
    # Using errorbar:
    plt.errorbar(res, dG_values, yerr=dG_errors, 
                 fmt='o-', 
                 color=cmap(i),
                 linewidth=2, 
                 markersize=5, 
                 capsize=3,      # Adds caps to error bars
                 elinewidth=1,   # Width of error bar lines
                 label=f'T = {temp}')

plt.xlabel('Residue ID', fontsize=14)
plt.ylabel('DG (kcal/mol)', fontsize=14)
plt.title('Residue ID vs Delta G with Jackknife Error (N=48)', fontsize=16)
plt.legend(fontsize=12)
plt.grid(True, alpha=0.3)
plt.ylim(-20, 30)
plt.tight_layout()
plt.savefig('residue_id_vs_deltaG_with_error.png', dpi=300, bbox_inches='tight')
print("Plot saved as 'residue_id_vs_deltaG_with_error.png'")

# Print residues with Delta G > 100
print("\n" + "="*50)
print("RESIDUES WITH DELTA G > 100")
print("="*50)
# ==========================================
# SAVE TO CSV
# ==========================================
# Converts the list of mean dG values to a numpy array and saves it.
# Rows = Temperatures, Columns = Residues.
# delimiter=',' makes it a CSV.
# No header or footer is added.
np.savetxt('glpg-dg.csv', np.array(dGhx_target_mean), delimiter=',')

for i, temp in enumerate(target_temps):
    dG_values = dGhx_target_mean[i]
    high_dG_mask = dG_values > 100
    high_dG_residues = res[high_dG_mask]
    
    print(f"\nTemperature T = {temp}:")
    if len(high_dG_residues) > 0:
        print(f"  {len(high_dG_residues)} residues with Delta G > 100:")
        print(f"    Residues: {', '.join(map(str, high_dG_residues))}")
    else:
        print("  No residues with Delta G > 100")
