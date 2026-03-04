import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
from scipy.interpolate import interp1d
import tables as  tb
import upside_engine as ue
from matplotlib.colors import BoundaryNorm
from matplotlib.colors import LogNorm
from matplotlib.colors import ListedColormap
from matplotlib.ticker import MaxNLocator

fbin1a = np.linspace(-14.,18.,161)
fbin1b = np.linspace(-18.,14.,161)
fc1a = (fbin1a[1:] + fbin1a[:-1])*0.5
fc1b = (fbin1b[1:] + fbin1b[:-1])*0.5
z = np.linspace(-10,20,121)

def compact_sigmoid(x, sharpness):
    y = x*sharpness;
    result = 0.25 * (y+2) * (y-1)**2
    result = np.where((y< 1), result, np.zeros_like(result))
    result = np.where((y>-1), result, np.ones_like (result))
    return result

lz = compact_sigmoid(z, 1.)
rz = 1-compact_sigmoid(z, 1.)

name = sys.argv[1]

th = 16.0


for name in sys.argv[1:]:

    with tb.open_file(name, 'r') as t:
        cb_energy       = t.root.cb_energy[:]
        icb_energy      = t.root.icb_energy[:]
        hb_energy       = t.root.hb_energy[:]
        seq             = t.root.names[:]
        cb_zmin         = t.root._v_attrs.cb_z_min
        cb_zmax         = t.root._v_attrs.cb_z_max
        hb_zmin         = t.root._v_attrs.hb_z_min
        hb_zmax         = t.root._v_attrs.hb_z_max
    
    xx = np.linspace(1,17,81)
    
    fig,axes = plt.subplots(nrows=4, ncols=5, figsize=(16,12))
    
    for i in range(4):
        for j in range(5):
            k = i*5+j
            s = seq[k]
    
            y11 = ue.clamped_spline_value(cb_energy[k,0], xx)
            y12 = ue.clamped_spline_value(cb_energy[k,1], xx)
            y13 = ue.clamped_spline_value(icb_energy[k],  xx)
        
            ax = axes[i,j]
    
            ax.plot(y11, '-')
            ax.plot(y12, '-')
            ax.plot(y13, '-')
            ax.grid()
            ax.set_title(s)
    
plt.show()
