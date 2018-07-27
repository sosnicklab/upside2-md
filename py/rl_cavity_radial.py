import sys, os
import cPickle as cp
import numpy as np
import tables as tb
from scipy.spatial.distance import pdist, cdist
sys.path.append(os.path.expanduser('~/upside/py/'))
from upside_config import chain_endpts

def write_cavity_radial(config_fn, ids, cavity_radii):
    assert ids.shape == cavity_radii.shape

    with tb.open_file(config_fn, 'a') as t:
        g = t.create_group(t.root.input.potential, 'cavity_radial')
        g._v_attrs.arguments = np.array(['pos'])

        t.create_array(g, 'id', obj=ids)
        t.create_array(g, 'radius', obj=cavity_radii)
        t.create_array(g, 'spring_constant', obj=np.ones(len(ids))*50.)

class gen_ub_placements:
    """Generate and apply cavity radius and unbound positions for complexes."""
    def __init__(self, init_fn, chain_breaks_fn, cavity_radius=0.):
        # load data
        with open(init_fn, 'rb') as f:
            self.xyz = cp.load(f)[...,0]

        with open(chain_breaks_fn) as f:
            chain_first_residue = np.array([int(i) for i in f.readline().split()])
            rl_chains = [int(i) for i in f.readline().split()]

        # slice out receptor and ligand positions
        n_res = self.xyz.shape[0]/3
        n_ch = rl_chains[0] + rl_chains[1]
        r_res_endpt = chain_endpts(n_res, chain_first_residue, rl_chains[0]-1)[1]
        l_res_endpt = chain_endpts(n_res, chain_first_residue, n_ch-1)[1]
        r_sele = np.arange(r_res_endpt*3)
        l_sele = np.arange(r_res_endpt*3, l_res_endpt*3)

        self.r_xyz = self.xyz[r_sele]
        self.l_xyz = self.xyz[l_sele]

        # center com
        self.r_xyz = self.r_xyz - self.r_xyz.mean(axis=0)
        self.l_xyz = self.l_xyz - self.l_xyz.mean(axis=0)

        # ToDo: cavity_radius should be max r dist to r_com + l_dist to l_com
        if cavity_radius:
            self.cavity_radius = cavity_radius
        else:
            # Base the radius on max distance in system
            self.cavity_radius = (pdist(self.xyz).max() + 10.)/2.

        # CA atom closest to com
        r_com_atom = np.argmin(cdist(self.r_xyz, np.zeros((1,3))))
        l_com_atom = np.argmin(cdist(self.l_xyz, np.zeros((1,3)))) + self.r_xyz.shape[0]
        r_com_res = r_com_atom / 3
        l_com_res = l_com_atom / 3
        self.r_com_atom = r_com_res*3 + 1
        self.l_com_atom = l_com_res*3 + 1

        # sample uniformly from cube
        samples = np.random.rand(2000, 3) - 0.5
        radial_dist = np.linalg.norm(samples, axis=1)
        # select those that fall in unit sphere
        samples = 2*samples[radial_dist < 0.5]

        # select those not nearer than 5 A from receptor
        samples = samples * (self.cavity_radius)
        rl_dist_min = np.zeros(len(samples))
        for i, smpl in enumerate(samples):
            rl_dist_min[i] = cdist(self.r_xyz, self.l_xyz + smpl).min()
        self.samples = samples[rl_dist_min > 5.]

        # set up displacements
        self.n_samples = self.samples.shape[0]
        self.xyz_frames = np.zeros((self.n_samples, self.xyz.shape[0], 3))
        self.xyz_frames[:, :self.r_xyz.shape[0], :] = self.r_xyz

        for i, smpl in enumerate(self.samples):
            self.xyz_frames[i, self.r_xyz.shape[0]:, :] = self.l_xyz + smpl

        # track sample usage
        self.counter = 0

    def apply(self, config_fn, make_unbound=True, scale_factor=1.0):
        with tb.open_file(config_fn, 'a') as t:
            if "output" in t.root:
                t.root.output._f_remove(recursive=True, force=True)
            if "cavity_radial" in t.root.input.potential:
                t.root.input.potential.cavity_radial._f_remove(recursive=True, force=True)

            # # testing: for visualizing samples
            # og = t.create_group(t.root, "output")
            # t.create_array(og, "pos", obj=self.xyz_frames[:,None,:,:])
            # t.create_array(og, "time", obj=np.arange(self.xyz_frames.shape[0]))

            if make_unbound:
                idx = self.counter % self.n_samples
                t.root.input.pos[:] = self.xyz_frames[idx,...,None]
        if make_unbound:
            self.counter += 1

        ids = np.array([self.r_com_atom, self.l_com_atom])
        cavity_radii = np.array([2.0, self.cavity_radius*scale_factor])
        write_cavity_radial(config_fn, ids, cavity_radii)