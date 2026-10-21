#!/usr/bin/env python
import sys
from string import ascii_lowercase
import tables
import numpy as np

H_bond=0.88
O_bond=1.24

# Chain labels
chain_labels = []
for c in ascii_lowercase:
    chain_labels.append(c.upper())
for c1 in ascii_lowercase:
    for c2 in ascii_lowercase:
        chain_labels.append(c1.upper() + c2.upper())

def vmag(x):# special version for systems
    assert x.shape[-2] == 3
    return np.sqrt(x[...,0,:]**2 + x[...,1,:]**2 + x[...,2,:]**2)

def vhat(x):  # special version for systems
    return x / vmag(x)[...,None,:]

def print_augmented_vtf(fname, sequence, traj, chain_first_residue):
    n_timestep, n_atom, three, n_system = traj.shape
    assert three == 3
    assert n_atom%3 ==0
    n_res = n_atom//3

    vtf = open(fname,'w')

    N  = traj[:,0::3].astype('f4')
    CA = traj[:,1::3].astype('f4')
    C  = traj[:,2::3].astype('f4')
    H  = N[:,1: ] - H_bond * vhat(vhat(C [:,:-1]-N[:,1: ]) + vhat(CA[:,1: ]-N[:,1: ]))
    O  = C[:,:-1] - O_bond * vhat(vhat(CA[:,:-1]-C[:,:-1]) + vhat(N [:,1: ]-C[:,:-1]))

    # H is placed from the previous C and O from the next N, so neither exists at a chain end;
    # proline (CPR is cis-proline) has no amide H
    chain_first = np.zeros(n_res, dtype=bool)
    chain_first[0] = True
    if chain_first_residue is not None:
        chain_first[chain_first_residue] = True
    chain_last = np.append(chain_first[1:], True)
    has_H = ~chain_first & ~np.isin(sequence, ['PRO','CPR'])
    has_O = ~chain_last
    chain_idx = np.cumsum(chain_first) - 1

    # write structure information
    atom_id = 0
    for ns in range(n_system):
        for nr in range(n_res):
            res = 'resid %i resname %s segid s%i chain %s' % (nr, sequence[nr], ns, chain_labels[chain_idx[nr]])

            vtf.write('atom %i name N  %s\n' % (atom_id+0, res))
            vtf.write('atom %i name CA %s\n' % (atom_id+1, res))
            vtf.write('atom %i name C  %s\n' % (atom_id+2, res))

            if not chain_first[nr]:
                vtf.write('bond %i:%i\n' % (prev_C,atom_id+0)) # prevC->N  bond
            vtf.write('bond %i:%i\n' % (atom_id+0,atom_id+1))  # N->CA bond
            vtf.write('bond %i:%i\n' % (atom_id+1,atom_id+2))  # CA->C bond

            H_id = atom_id+3
            O_id = atom_id+3+has_H[nr]
            if has_H[nr]:
                vtf.write('atom %i name H  %s\n' % (H_id, res))
            if has_O[nr]:
                vtf.write('atom %i name O  %s\n' % (O_id, res))
            if has_H[nr]:
                vtf.write('bond %i:%i\n' % (atom_id+0, H_id))  # N->H bond
            if has_O[nr]:
                vtf.write('bond %i:%i\n' % (atom_id+2, O_id))  # C->O bond
            consumed = 3+has_H[nr]+has_O[nr]

            prev_C = atom_id+2
            atom_id += consumed

    for f in range(len(traj)):
        vtf.write("\ntimestep ordered\n")
        for ns in range(n_system):
            for nr in range(n_res):
                vtf.write("%.3f %.3f %.3f\n" % (N [f,nr,0,ns], N [f,nr,1,ns], N [f,nr,2,ns]))
                vtf.write("%.3f %.3f %.3f\n" % (CA[f,nr,0,ns], CA[f,nr,1,ns], CA[f,nr,2,ns]))
                vtf.write("%.3f %.3f %.3f\n" % (C [f,nr,0,ns], C [f,nr,1,ns], C [f,nr,2,ns]))
                if has_H[nr]:
                    vtf.write("%.3f %.3f %.3f\n" % (H[f,nr-1,0,ns], H[f,nr-1,1,ns], H[f,nr-1,2,ns]))
                if has_O[nr]:
                    vtf.write("%.3f %.3f %.3f\n" % (O[f,nr,0,ns], O[f,nr,1,ns], O[f,nr,2,ns]))

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input_h5', help='Input simulation file')
    parser.add_argument('output_vtf', help='Output trajectory file')
    parser.add_argument('--stride', type=int, default=1, help='Stride for reading file')
    parser.add_argument('--start', type=int, default=0, help='Initial frame to extract VTF')
    parser.add_argument('--top_h5', type=str, default=None, help='Input top file')
    args = parser.parse_args()

    top_file = args.input_h5
    if args.top_h5:
        top_file = args.top_h5

    with tables.open_file(top_file) as t:
        seq = np.char.decode(t.root.input.sequence[:], encoding='UTF-8')
        if 'chain_break' in t.root.input:
            chain_first_residue = t.root.input.chain_break.chain_first_residue[:]
        else:
            chain_first_residue = None

    with tables.open_file(args.input_h5) as t:
        output_paths = []
        i = 0
        while 'output_previous_%i'%i in t.root:
            output_paths.append('/output_previous_%i'%i)
            i+=1
        if 'output' in t.root:  # 'output' is the *last* produced output
            output_paths.append('/output')

        total_frames_produced = 0
        pos = []
        stride = args.stride
        for g_no, opath in enumerate(output_paths):
            g = t.get_node(t.get_node(opath))
            # the first frame of each pos after the first is the same as the last frame before restart;
            # take the frames args.start + k*stride of the concatenated trajectory
            first = 1 if g_no else 0
            start_frame = max(args.start, total_frames_produced)
            start_frame += (args.start - start_frame) % stride
            pos.append(g.pos[start_frame-total_frames_produced+first::stride].transpose((0,2,3,1)))
            total_frames_produced += g.pos.shape[0]-first
            print (opath, total_frames_produced, 'cumulative frames found')

    pos = np.concatenate(pos, axis=0)
    print (pos.shape[0], 'frames for output', pos.shape[0]/30., 'seconds at 30 frames/second video')
    print_augmented_vtf(args.output_vtf, seq, pos, chain_first_residue)

if __name__ == '__main__':
    main()
