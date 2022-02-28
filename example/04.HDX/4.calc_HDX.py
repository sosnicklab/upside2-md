import numpy as np
import matplotlib.pyplot as plt
import pymbar # for MBAR analysis
from pymbar import timeseries # for timeseries analysis

pdb_id      = 'EHEE_rd2_0005'
sim_id      = 'REMD'
start_frame = 100

work_dir         = './'
n_rep            = 16     # replica number

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
            print(k, r)
            dG[r] = 1000.
        else:
            dG[r] = 0.001987*tt*np.log((mean_pf/(1.-mean_pf)))
    dGhx_T.append(dG)
dGhx_T = np.array(dGhx_T)

for i,r in enumerate(res):
    plt.plot(T, dGhx_T[:,i])

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
# plot
#============================================

# [den] vs DG
fig = plt.figure()
for i,r in enumerate(res):
    plt.plot(den, dGhx_D[:,i])
plt.xlabel('[den]')
plt.ylabel('DG (kcal/mol)')

# DG of residues
fig = plt.figure()
plt.bar(res, dGhx_D[0,:]) 
plt.xlabel('seq')
plt.ylabel('DG (kcal/mol)')
plt.show()
