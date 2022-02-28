import numpy as np
import matplotlib.pyplot as plt
import pymbar # for MBAR analysis
from pymbar import timeseries # for timeseries analysis

pdb_id      = 'EHEE_rd2_0005'
sim_id      = 'REMD'
start_frame = 50

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
for i in range(n_rep):
    Pot.append(  np.load('{}/{}_{}_{}_Energy.npy'.format(result_dir, pdb_id, sim_id, i))[:,0] ) 
    Rg.append(   np.load('{}/{}_{}_{}_Rg.npy'    .format(result_dir, pdb_id, sim_id, i)) ) 
    Hb.append(   np.load('{}/{}_{}_{}_Hbond.npy' .format(result_dir, pdb_id, sim_id, i)) ) 
    Rmsd.append( np.load('{}/{}_{}_{}_Rmsd.npy'  .format(result_dir, pdb_id, sim_id, i)) ) 
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
# reweighting
#============================================


# melting cure

Rg_mean = []
Hb_mean = []
Rg_reweight = []
Hb_reweight = []

Rg_flatten = Rg[:,start_frame:].flatten()
Hb_flatten = Hb[:,start_frame:].flatten()

for k in range(T.size):
    t = T[k]
    u_n = (cE0/(t*kB)).flatten()
    log_w1 = mbar0._computeUnnormalizedLogWeights(u_n)
    w1 = np.exp(log_w1)
    w1 /= np.sum(w1)

    Rg_mean.append( np.mean( Rg[k, start_frame:]) )
    Hb_mean.append( np.mean( Hb[k, start_frame:]) )
    Rg_reweight.append( np.sum( Rg_flatten * w1) )
    Hb_reweight.append( np.sum( Hb_flatten * w1) )

# free energy along the #Hbond
fig1 = plt.figure()
plt.plot(T, Rg_mean, label='direct average' )
plt.plot(T, Rg_reweight, label='MBAR reweighting' )
plt.legend()
plt.xlabel('T')
plt.ylabel('Rg')

dG_hbond = []
for k in range(T.size):
    t = T[k]
    u_n = (cE0/(t*kB)).flatten()
    log_w1 = mbar0._computeUnnormalizedLogWeights(u_n)
    w1 = np.exp(log_w1)
    w1 /= np.sum(w1)

    a,b = np.histogram(Hb_flatten, 20, weights=w1)
    bc = (b[1:]+b[:-1])*0.5

    g = -0.593*np.log(a)*t
    g -= np.min(g)
    dG_hbond.append(g)

fig2 = plt.figure()
selected_replica = 8
plt.plot(bc, dG_hbond[selected_replica] )
plt.title('T={:.2f}'.format(T[selected_replica]))
plt.xlabel('#H-bond')
plt.ylabel('Free Energy (kcal/mol)')

plt.show()
