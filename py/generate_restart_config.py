#!/usr/bin/env python
'''
Generate a new config file (.up or .h5) for restart simulation.
Run it with --restart-using-momentum when the source run recorded momentum (--record-momentum).
'''

__author__  = 'Xiangda Peng'
__version__ = '2020-05-11'

import sys,os
import numpy as np
import tables as tb
import h5py as h5


def generate_config_for_restart(up_fname, up_fname_for_next_round, frame=-1):

    # copy the root/input group from up_fname to the new up
    with h5.File(up_fname, 'r') as in_H5:
        in_input = in_H5['/input']
        with h5.File(up_fname_for_next_round, 'a') as out_H5:
            if not ('input' in out_H5.keys()):
                out_root = out_H5['/']
                in_H5.copy(in_input, out_root)

    # write the last frame to a new up file
    with tb.File(up_fname, 'r') as in_H5:
        with tb.File(up_fname_for_next_round, 'a') as out_H5:
            NF = in_H5.root.output.pos.shape[0]
            assert NF > frame

            out_H5.root.input.pos[:,:,0] = in_H5.root.output.pos[frame, 0]
            # the momentum of the frame and the round it was recorded at, which --restart-using-momentum
            # reads to resume the momentum and every random counter; a run without --record-momentum has none
            if 'mom' in out_H5.root.input:
                out_H5.root.input.mom._f_remove()
            if 'mom' in in_H5.root.output:
                out_H5.create_array(out_H5.root.input, 'mom', obj=in_H5.root.output.mom[frame, 0][:,:,None])
                out_H5.root.input.mom._v_attrs.round = in_H5.root.output.round[frame]
                if 'replica_index' in in_H5.root.output:
                    # the replica the system held at the frame, from which replica exchange continues
                    out_H5.root.input.mom._v_attrs.replica_index = in_H5.root.output.replica_index[frame, 0]
            if 'box' in in_H5.root.output:
                # the barostat's box at the frame, in every node that carries one
                box = [float(v) for v in in_H5.root.output.box[frame]]
                for node in out_H5.root.input.potential:
                    if all(k in node._v_attrs for k in ('x_len', 'y_len', 'z_len')):
                        node._v_attrs.x_len, node._v_attrs.y_len, node._v_attrs.z_len = box
            if 'curvature_center' in in_H5.root.output:
                # --curvature-changer-interval moves the centre of curvature; continue from the frame's
                out_H5.root.input.potential.Const3D_curvature_center.value[:] = in_H5.root.output.curvature_center[frame]
            if 'tip_pos' in in_H5.root.output:
                # the tip is at start_pos + velocities*(time_initial + time_step*steps): resume its step count
                out_H5.root.input.potential.MovingConst3D._v_attrs.step_initial = in_H5.root.output.tip_step[frame, 0]

def main():
    if len(sys.argv) not in (3, 4):
        raise SystemExit('usage: generate_restart_config.py CONFIG CONFIG_FOR_NEXT_ROUND [FRAME]')
    up_fname                = sys.argv[1]
    up_fname_for_next_round = sys.argv[2]
    frame = int(sys.argv[3]) if len(sys.argv)==4 else -1
    generate_config_for_restart(up_fname, up_fname_for_next_round, frame)

if __name__ == '__main__':
    main()
