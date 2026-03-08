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

name = sys.argv[1]

for name in sys.argv[1:]:

    with tb.open_file(name, 'r') as t:
        cb_energy       = t.root.cb_energy[:]
        icb_energy      = t.root.icb_energy[:]
        hb_energy       = t.root.hb_energy[:]
        ihb_energy      = t.root.ihb_energy[:]
        seq             = t.root.names[:]
        cb_zmin         = t.root._v_attrs.cb_z_min
        cb_zmax         = t.root._v_attrs.cb_z_max
        hb_zmin         = t.root._v_attrs.hb_z_min
        hb_zmax         = t.root._v_attrs.hb_z_max
    
    xx = np.linspace(1,17,81)
    
    fig,axes = plt.subplots(nrows=1, ncols=2, figsize=(16,12))

    for i in [0,1]:

        ax = axes[i]
        
        y11 = ue.clamped_spline_value(hb_energy[i,0], xx)
        y12 = ue.clamped_spline_value(hb_energy[i,1], xx)
        y13 = ue.clamped_spline_value(ihb_energy[i,0], xx)
        y14 = ue.clamped_spline_value(ihb_energy[i,1], xx)
            
        ax.plot(y11, 'r-')
        ax.plot(y12, 'b-')
        ax.plot(y13, 'r--')
        ax.plot(y14, 'b--')
        ax.grid()
    
plt.show()
