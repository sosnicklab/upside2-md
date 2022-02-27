#!/usr/bin/env python
'''
Generate a new config file (.up or .h5) for restart simulation
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
            if 'tip_pos' in in_H5.root.output:
                use_AFM_pulling = True
                tip_pos         = in_H5.root.output.tip_pos[frame]
                time_estimate   = in_H5.root.output.time_estimate[frame][0]
                out_H5.root.input.potential.AFM.starting_tip_pos[:] = tip_pos
                #out_H5.root.input.potential.AFM.pulling_vel._v_attrs.time_initial = time_estimate
            if 'rotation_tip_angles' in in_H5.root.output:
                use_rotation = True
                rotation_tip_angles    = in_H5.root.output.rotation_tip_angles[frame]
                rotation_time_estimate = in_H5.root.output.rotation_time_estimate[frame][0]
                out_H5.root.input.potential.rotation.starting_tip_angle[:] = rotation_tip_angles
                #out_H5.root.input.potential.rotation.angular_vel._v_attrs.time_initial = rotation_time_estimate

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
