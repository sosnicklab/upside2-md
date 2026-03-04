import cPickle as cp
import numpy as np
import sys

data     = cp.load(open('{}.initial.pkl'.format(sys.argv[1])))
xyz      = data[:,:,0]
n_res    = int(xyz.shape[0]/3.)
ndx_ca   = np.arange(n_res)*3+1
ndx_tm   = np.where((xyz[ndx_ca,2]<15.3)*(xyz[ndx_ca,2]>-15.3))[0]
ndx_sres = np.random.randint(ndx_tm.size, size=20)

print "residue spring_const"
for n in ndx_sres:
    print ndx_tm[n], 4.0

