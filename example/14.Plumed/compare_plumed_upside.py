#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import mdtraj_upside as mu
import os
import matplotlib.pyplot as plt


# ## compare coordinates

# In[2]:


traj = mu.load_upside_traj("outputs/simple_test/chig.run.up")


# In[3]:


def quick_xyz_parser(file, natoms):
    coords_all = []
    i=0
    with open(file, 'r') as f:
        for line in f.readlines():
            if line.rstrip("\n") == str(natoms):
                i+=1
                try:
                    coords
                except NameError:
                    coords = []
                else:
                    coords_all.append(coords)
                    coords = []
            elif line[0] != 'X':
                pass
            else:
                X, Y, Z = line.lstrip("X ").rstrip("\n").split(" ")

                coords.append([float(X), float(Y), float(Z)])
        coords_all.append(coords)
    return coords_all


# In[5]:


coords_all = np.array(quick_xyz_parser('outputs/simple_test/chig.xyz.plumed', natoms=30))
coords_all.shape


# In[6]:


core_ids = traj.topology.select("name N or name CA or name C")


# In[7]:


traj.xyz[:, core_ids].shape


# In[10]:


if not os.path.exists("results"):
    os.mkdir("results")

plt.plot(traj.time, np.sum((coords_all - traj.xyz[:, core_ids])**2, axis=(1,2)))
plt.xlabel("Time (in Upside unit)")
plt.ylabel("RMS difference in coordinates")
plt.savefig("results/compare_pos_upside_plumed.png", dpi=200, bbox_inches='tight')
plt.show()


# ## compare Rg

# In[11]:


def calc_rg(xyz):
    return np.sqrt(np.sum((xyz - xyz.mean(axis=(0,1)))**2)/ len(xyz[0, :]))


# In[12]:


rgs_upside = [calc_rg(t.xyz[:, core_ids]) for t in traj]
rgs_upside


# In[14]:


rgs_plumed = np.loadtxt("outputs/simple_test/chig.rg.plumed")[:, 1]


# In[16]:


plt.plot(traj.time, rgs_upside, '-', label='mdtraj calculation')
plt.plot(traj.time, rgs_plumed, '--', label='plumed calculation')
plt.legend()
plt.xlabel("Time (in Upside unit)")
plt.ylabel("Rg (nm)")
plt.savefig("results/compare_Rg_upside_plumed.png", dpi=200, bbox_inches='tight')
plt.show()


# In[ ]:




