#!/usr/bin/env python
# coding: utf-8




import numpy as np
import mdtraj_upside as mu
import mdtraj as md
import os
import matplotlib.pyplot as plt

POS = False
RG = False
RMSD = True

# ## compare coordinates




#traj = mu.load_upside_traj("outputs/simple_test/chig.run.up")
traj = mu.load_upside_traj("outputs/metad/chig.run.up")





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


if not os.path.exists("results"):
    os.mkdir("results")



if POS:
    coords_all = np.array(quick_xyz_parser('outputs/simple_test/chig.xyz.plumed', natoms=30))
    core_ids = traj.topology.select("name N or name CA or name C")

    plt.plot(traj.time, np.sum((coords_all - traj.xyz[:, core_ids])**2, axis=(1,2)))
    plt.xlabel("Time (in Upside unit)")
    plt.ylabel("RMS difference in coordinates")
    plt.savefig("results/compare_pos_upside_plumed.png", dpi=200, bbox_inches='tight')
    plt.show()


# ## compare Rg

def calc_rg(xyz):
    return np.sqrt(np.sum((xyz - xyz.mean(axis=(0,1)))**2)/ len(xyz[0, :]))

if RG:
    rgs_upside = [calc_rg(t.xyz[:, core_ids]) for t in traj]
    rgs_plumed = np.loadtxt("outputs/simple_test/chig.rg.plumed")[:, 1]
    plt.plot(traj.time, rgs_upside, '-', label='mdtraj calculation')
    plt.plot(traj.time, rgs_plumed, '--', label='plumed calculation')
    plt.legend()
    plt.xlabel("Time (in Upside unit)")
    plt.ylabel("Rg (nm)")
    plt.savefig("results/compare_Rg_upside_plumed.png", dpi=200, bbox_inches='tight')
    plt.show()


# ## Compare RMSD 
if RMSD:
    sele = traj.top.select("name CA")
    rmsd_upside = md.rmsd(traj, traj, atom_indices=sele)
    rmsd_plumed = np.loadtxt("outputs/metad/chig.rmsd.plumed")[:,1]
    plt.plot(rmsd_upside, '-', label='mdtraj calculation')
    plt.plot(rmsd_plumed, '--', label='plumed calculation')
    plt.legend()
    plt.xlabel("Time (in Upside unit)")
    plt.ylabel("Rmsd (nm)")
    plt.savefig("results/compare_Rmsd_upside_plumed.png", dpi=200, bbox_inches='tight')
    plt.show()






