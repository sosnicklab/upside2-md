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

fbin1a = np.linspace(-14.,20.,171)
fbin1b = np.linspace(-20.,14.,171)
fc1a = (fbin1a[1:] + fbin1a[:-1])*0.5
fc1b = (fbin1b[1:] + fbin1b[:-1])*0.5
z = np.linspace(-30,30,61)

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

with tb.open_file(name, 'r') as t:
    cb_energy_left  = t.root.cb_energy_left[:]
    cb_energy_right = t.root.cb_energy_right[:]
    hb_energy_left  = t.root.uhb_energy_left[:]
    hb_energy_right = t.root.uhb_energy_right[:]
    burial_nodes    = t.root.burial_nodes[:]
    z_min_left      = t.root._v_attrs.z_min_left
    z_max_left      = t.root._v_attrs.z_max_left
    z_min_right     = t.root._v_attrs.z_min_right
    z_max_right     = t.root._v_attrs.z_max_right
    seq             = t.root.names[:]

clz = z-(z_min_left-th)
crz = z-(z_min_right+th)

for i,s in enumerate(seq):

    y11 = ue.clamped_spline_value(cb_energy_left[i,0],  clz)
    y12 = ue.clamped_spline_value(cb_energy_left[i,1],  clz)
    y13 = ue.clamped_spline_value(cb_energy_left[i,2],  clz)
    y14 = ue.clamped_spline_value(cb_energy_left[i,3],  clz)
    y15 = ue.clamped_spline_value(cb_energy_left[i,4],  clz)
    
    y21 = ue.clamped_spline_value(cb_energy_right[i,0], crz)
    y22 = ue.clamped_spline_value(cb_energy_right[i,1], crz)
    y23 = ue.clamped_spline_value(cb_energy_right[i,2], crz)
    y24 = ue.clamped_spline_value(cb_energy_right[i,3], crz)
    y25 = ue.clamped_spline_value(cb_energy_right[i,4], crz)

    y1 = lz*y11 + rz*y21
    y2 = lz*y12 + rz*y22
    y3 = lz*y13 + rz*y23
    y4 = lz*y14 + rz*y24
    y5 = lz*y15 + rz*y25

    fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(10,3))
    ax.plot(z, y1, '-')
    ax.plot(z, y2, '-')
    ax.plot(z, y3, '-')
    ax.plot(z, y4, '-')
    ax.plot(z, y5, '-')
    ax.set_title(s)
    
    #fig,axes = plt.subplots(nrows=1, ncols=2, figsize=(10,3))
    #ax = axes[0]
    #ax.plot(fc1a, y11, '-')
    #ax.plot(fc1a, y12, '-')
    #ax.plot(fc1a, y13, '-')
    #ax.plot(fc1a, y14, '-')
    #ax.plot(fc1a, y15, '-')
    #ax.set_title(s)
    #
    #ax = axes[1]
    #ax.plot(fc1b, y21, '-')
    #ax.plot(fc1b, y22, '-')
    #ax.plot(fc1b, y23, '-')
    #ax.plot(fc1b, y24, '-')
    #ax.plot(fc1b, y25, '-')

#    plt.savefig('{}-MZ.png'.format(res))
plt.show()
