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
            if 'tip_pos' in in_H5.root.output:
                # the tip is at start_pos + velocities*time and its clock restarts at time_initial,
                # so keep the tip's start and resume its clock
                moving        = out_H5.root.input.potential.MovingConst3D
                time_estimate = in_H5.root.output.time_estimate[frame][0]
                start_pos     = in_H5.root.output.tip_pos[frame] - moving.velocities[:]*time_estimate
                if 'start_pos' in moving:
                    moving.start_pos._f_remove()
                out_H5.create_array(moving, 'start_pos', obj=start_pos)
                moving._v_attrs.initialized_by_coord = 0
                moving._v_attrs.time_initial         = time_estimate

def main():
    up_fname                = sys.argv[1]
    up_fname_for_next_round = sys.argv[2]
    if len(sys.argv)==3:
        generate_config_for_restart(up_fname, up_fname_for_next_round)
    elif len(sys.argv)==4:
        frame = int(sys.argv[3])
        generate_config_for_restart(up_fname, up_fname_for_next_round, frame)

if __name__ == '__main__':
    main()
