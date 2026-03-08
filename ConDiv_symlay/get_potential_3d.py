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

bin1a = np.linspace(-13., 17., 16)
bin1b = np.linspace(-17., 13., 16)
c1a = (bin1a[1:] + bin1a[:-1])*0.5
c1b = (bin1b[1:] + bin1b[:-1])*0.5

dx = c1a[1] - c1a[0]
zmin_left  = c1a[0] - dx
zmin_right = c1b[0] - dx

fbin1a = np.linspace(-14.,14.,141)
fbin1b = np.linspace(-14.,14.,141)
fc1a = (fbin1a[1:] + fbin1a[:-1])*0.5
fc1b = (fbin1b[1:] + fbin1b[:-1])*0.5

f = '/home/pengxd/upside-ff2.0/parameters/membrane_potential/UC_mempot_ref-surf_exposed_thickness38.0_unit-RT.h5'
with tb.open_file(f) as t:
    z_min = t.root.uhb_energy._v_attrs.z_min
    z_max = t.root.uhb_energy._v_attrs.z_max
    unh   = t.root.uhb_energy[0]
    uco   = t.root.uhb_energy[1]
z = np.linspace(z_min, z_max, 501)

th = 19
z1 = z+th
z2 = z-th

f11 = interp1d(z1, unh)
f12 = interp1d(z1, uco)
f21 = interp1d(z2, unh)
f22 = interp1d(z2, uco)

py11 = f11(c1a)
py12 = f12(c1a)
py21 = f21(c1b)
py22 = f22(c1b)

p11 = ue.clamped_spline_solve(py11)
p12 = ue.clamped_spline_solve(py12)
p21 = ue.clamped_spline_solve(py21)
p22 = ue.clamped_spline_solve(py22)

n_node = p11.size

hb_param_left  = np.zeros((2, n_node))
hb_param_right = np.zeros((2, n_node))

hb_param_left[0]  = p11
hb_param_left[1]  = p12
hb_param_right[0] = p21
hb_param_right[1] = p22

y11 = ue.clamped_spline_value(p11, (fc1a-zmin_left)*0.5)
y12 = ue.clamped_spline_value(p12, (fc1a-zmin_left)*0.5)
y21 = ue.clamped_spline_value(p21, (fc1b-zmin_right)*0.5)
y22 = ue.clamped_spline_value(p22, (fc1b-zmin_right)*0.5)

fig,axes = plt.subplots(nrows=1, ncols=2, figsize=(10,3))
ax = axes[0]
ax.plot(c1a, py11)
ax.plot(c1a, py12)
ax.plot(fc1a, y11, '--')
ax.plot(fc1a, y12, '--')

ax = axes[1]
ax.plot(c1b, py21)
ax.plot(c1b, py22)
ax.plot(fc1b, y21, '--')
ax.plot(fc1b, y22, '--')
plt.show()

seq = np.loadtxt('seq', dtype=str)

data1    = np.loadtxt('left')
ndx_1    = np.where((data1[:,0]>-16)*(data1[:,0]<30))[0]
kernel01 = st.gaussian_kde(data1[ndx_1,0])
f01     = kernel01(c1a)
f01    /= np.sum(f01)

data2    = np.loadtxt('right')
ndx_1    = np.where((data2[:,0]>-30)*(data2[:,0]<16))[0]
kernel02 = st.gaussian_kde(data2[ndx_1,0])
f02     = kernel02(c1b)
f02    /= np.sum(f02)

def compact_sigmoid(x, sharpness):
    y = x*sharpness;
    result = 0.25 * (y+2) * (y-1)**2
    result = np.where((y< 1), result, np.zeros_like(result))
    result = np.where((y>-1), result, np.ones_like (result))
    return result

bl_node = [2.0, 4.0]
n_bl = len(bl_node) + 1
t
param_left  = np.zeros((20, n_bl, n_node))
param_right = np.zeros((20, n_bl, n_node))

for i,res in enumerate(seq[:]):

    data1    = np.loadtxt('{}.left'.format(res))
    ndx_1    = np.where((data1[:,0]>-16)*(data1[:,0]<30)*(data1[:,1]<2.0))[0]
    ndx_2    = np.where((data1[:,0]>-16)*(data1[:,0]<30)*(data1[:,1]>=2.0)*(data1[:,1]<4.0))[0]
    ndx_3    = np.where((data1[:,0]>-16)*(data1[:,0]<30)*(data1[:,1]>=4.0))[0]

    kernel11 = st.gaussian_kde(data1[ndx_1,0])
    kernel12 = st.gaussian_kde(data1[ndx_2,0])
    kernel13 = st.gaussian_kde(data1[ndx_3,0])
    f11      = kernel11(c1a)
    f12      = kernel12(c1a)
    f13      = kernel13(c1a)
    f11      /= np.sum(f11)/ndx_1.size
    f12      /= np.sum(f12)/ndx_2.size
    f13      /= np.sum(f13)/ndx_3.size

    data2    = np.loadtxt('{}.right'.format(res))
    ndx_1    = np.where((data2[:,0]>-30)*(data2[:,0]<16)*(data2[:,1]<2.0))[0]
    ndx_2    = np.where((data2[:,0]>-30)*(data2[:,0]<16)*(data2[:,1]>=2.0)*(data2[:,1]<4.0))[0]
    ndx_3    = np.where((data2[:,0]>-30)*(data2[:,0]<16)*(data2[:,1]>=4.0))[0]

    kernel21 = st.gaussian_kde(data2[ndx_1,0])
    kernel22 = st.gaussian_kde(data2[ndx_2,0])
    kernel23 = st.gaussian_kde(data2[ndx_3,0])
    f21      = kernel21(c1b)
    f22      = kernel22(c1b)
    f23      = kernel23(c1b)
    f21      /= np.sum(f21)/ndx_1.size
    f22      /= np.sum(f22)/ndx_2.size
    f23      /= np.sum(f23)/ndx_3.size

    g11 = f11/f01
    g12 = f12/f01
    g13 = f13/f01

    g21 = f21/f02
    g22 = f22/f02
    g23 = f23/f02

    g11 = -np.log(g11)
    g12 = -np.log(g12)
    g13 = -np.log(g13)
    g21 = -np.log(g21)
    g22 = -np.log(g22)
    g23 = -np.log(g23)

    ap1 = np.where(c1a==14)[0][0]
    ap2 = np.where(c1b==-14)[0][0]
    g21 += g11[ap1]-g21[ap2]
    g22 += g12[ap1]-g22[ap2]
    g23 += g13[ap1]-g23[ap2]

    minv = np.min([g11, g12, g13, g21, g22, g23])
    maxv = np.max([g11, g12, g13, g21, g22, g23])

    g11 -= minv
    g12 -= minv
    g13 -= minv
    g21 -= minv
    g22 -= minv
    g23 -= minv

    minv = np.min([g11, g12, g13, g21, g22, g23])
    maxv = np.max([g11, g12, g13, g21, g22, g23])

    p11 = ue.clamped_spline_solve(g11)
    p12 = ue.clamped_spline_solve(g12)
    p13 = ue.clamped_spline_solve(g13)

    p21 = ue.clamped_spline_solve(g21)
    p22 = ue.clamped_spline_solve(g22)
    p23 = ue.clamped_spline_solve(g23)

    param_left[i,0]  = p11
    param_left[i,1]  = p12
    param_left[i,2]  = p13

    param_right[i,0] = p21
    param_right[i,1] = p22
    param_right[i,2] = p23

    y11 = ue.clamped_spline_value(p11, (fc1a-zmin_left)*0.5)
    y12 = ue.clamped_spline_value(p12, (fc1a-zmin_left)*0.5)
    y13 = ue.clamped_spline_value(p13, (fc1a-zmin_left)*0.5)

    y21 = ue.clamped_spline_value(p21, (fc1b-zmin_right)*0.5)
    y22 = ue.clamped_spline_value(p22, (fc1b-zmin_right)*0.5)
    y23 = ue.clamped_spline_value(p23, (fc1b-zmin_right)*0.5)

    fig,axes = plt.subplots(nrows=1, ncols=2, figsize=(10,3))
    ax = axes[0]
    ax.set_title(res)
    ax.plot(c1a, g11)
    ax.plot(c1a, g12)
    ax.plot(c1a, g13)

    ax.plot(fc1a, y11, '--')
    ax.plot(fc1a, y12, '--')
    ax.plot(fc1a, y13, '--')

    ax.set_title(res)
    ax.plot([14,14], [minv, maxv], 'k--')
    ax.plot([0,0], [minv, maxv], 'k--')
    ax = axes[1]
    ax.plot(c1b, g21)
    ax.plot(c1b, g22)
    ax.plot(c1b, g23)

    ax.plot(fc1b, y21, '--')
    ax.plot(fc1b, y22, '--')
    ax.plot(fc1b, y23, '--')

    ax.plot([-14,-14], [minv, maxv], 'k--')
    ax.plot([0,0], [minv, maxv], 'k--')
    plt.tight_layout()
    plt.savefig('{}-MZ.png'.format(res))

with tb.open_file('membrane.h5', 'a') as t:
    root = t.root
    t.create_array(root, 'cb_energy_left',   obj=param_left)
    t.create_array(root, 'cb_energy_right',  obj=param_right)
    t.create_array(root, 'uhb_energy_left',  obj=hb_param_left)
    t.create_array(root, 'uhb_energy_right', obj=hb_param_right)
    t.create_array(root, 'names',            obj=seq)
    t.create_array(root, 'burial_nodes',     obj=np.array(bl_node))

    t.root._v_attrs.z_min_left = zmin_left
    t.root._v_attrs.z_max_left = c1a[-1]
    t.root._v_attrs.z_min_right = zmin_right
    t.root._v_attrs.z_max_right = c1b[-1]
    print t.root._v_attrs.z_min_left, t.root._v_attrs.z_max_left, t.root._v_attrs.z_min_right, t.root._v_attrs.z_max_right, t.root.cb_energy_left.shape

