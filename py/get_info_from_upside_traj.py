import sys, os
upside_path = os.environ['UPSIDE_HOME']
upside_utils_dir = os.path.expanduser(upside_path+"/py")
sys.path.insert(0, upside_utils_dir)

import mdtraj as md
import mdtraj_upside as mu
import upside_engine as ue
import numpy as np
import tables as tb

def _output_groups(t):
    i=0
    while 'output_previous_%i'%i in t.root:
        yield t.get_node('/output_previous_%i'%i)
        i += 1
    if 'output' in t.root:
        yield t.get_node('/output')
        i += 1

def _output_groups_reverse(t, n):
    if 'output' in t.root:
        yield t.get_node('/output')

    for i in range(n-2,-1,-1):
        if 'output_previous_%i'%i in t.root:
            yield t.get_node('/output_previous_%i'%i)

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input_h5', help='Input simulation file')
    parser.add_argument('output_base', help='Output name')
    parser.add_argument('--top', type=str, default=None, help='Input top/ref file (.up or .h5 file)')
    parser.add_argument('--stride', type=int, default=1, help='(default 1) Stride for reading file')
    parser.add_argument('--start', type=int, default=0, help='(default 0) the first frame')
    parser.add_argument('--last', type=int, default=None, help='(default none) the last frame')
    args = parser.parse_args()

    first  = args.start
    stride = args.stride
    last   = args.last

    # use mdtraj_upside to load the trajs
    if (args.top):
        trj = mu.load_upside_traj(args.input_h5, top=args.top, add_atoms=False)
        re  = mu.load_upside_ref(args.top, add_atoms=False)
    else:
        trj = mu.load_upside_traj(args.input_h5, add_atoms=False)
        re = trj[0]

    if not args.last:
        last = trj.n_frames

    sele   = trj.top.select("name CA")
    Rmsd = 10.*md.rmsd(trj, re, atom_indices=sele)[first:last:stride]
    Rg   = 10.*md.compute_rg(trj)[first:last:stride]

    # use the "tables" to open the h5 file directly
    exchange = False
    with tb.open_file(args.input_h5) as t:
        T   = t.root.output.temperature[0,0]
     
        for g_no, g in enumerate(_output_groups(t)):
            # take into account that the first frame of each pos is the same as the last frame before restart
            # attempt to land on the stride
            if g_no == 0:
                Pot = g.potential[:]
                Hb  = np.sum(g.hbond[:], axis=1)
            else:
                Pot = np.concatenate([ Pot, g.potential[1:] ]) 
                Hb  = np.concatenate([ Hb,  np.sum(g.hbond[1:], axis=1) ])
        Pot = Pot[first:last:stride]
        Hb  = Hb[first:last:stride]

    np.save('{}_Energy.npy'.format(args.output_base), Pot )
    np.save('{}_Hbond.npy'.format(args.output_base), Hb )
    np.save('{}_Rmsd.npy'.format(args.output_base), Rmsd )
    np.save('{}_Rg.npy'.format(args.output_base), Rg )
    np.save('{}_T.npy'.format(args.output_base), T )

if __name__ == '__main__':
    main()
