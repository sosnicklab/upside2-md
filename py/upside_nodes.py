import numpy as np
import tables as tb
import sys,os

deg=np.deg2rad(1)
default_filter = tb.Filters(complib='zlib', complevel=5, fletcher32=True)

def bstring(string):
    return bytes(string, encoding="ascii")

def create_array(h5, grp, nm, obj=None):
    return h5.create_earray(grp, nm, obj=obj, filters=default_filter)

def DistanceCoord(h5, nodename, input_node1, input_node2, ids ):
    potential = h5.root.input.potential
    base = 'Distance3D'
    if nodename is None:
        name = base
    else:
        name = '{}_{}'.format(base, nodename)
    dgrp = h5.create_group(potential, name)
    dgrp._v_attrs.arguments = np.array([bstring(input_node1), bstring(input_node2)])
    create_array(h5, dgrp, 'id', obj=ids)
    return name

def Distance2DCoord(h5, nodename, input_node1, input_node2, dim1, dim2, ids ):
    potential = h5.root.input.potential
    base = 'Distance2D'
    if nodename is None:
        name = base
    else:
        name = '{}_{}'.format(base, nodename)
    dgrp = h5.create_group(potential, name)
    dgrp._v_attrs.arguments = np.array([bstring(input_node1), bstring(input_node2)])
    dgrp._v_attrs.dim1 = dim1
    dgrp._v_attrs.dim2 = dim2
    create_array(h5, dgrp, 'id', obj=ids)
    return name

def Distance1DCoord(h5, nodename, input_node1, input_node2, dim1, ids ):
    potential = h5.root.input.potential
    base = 'Distance1D'
    if nodename is None:
        name = base
    else:
        name = '{}_{}'.format(base, nodename)
    dgrp = h5.create_group(potential, name)
    dgrp._v_attrs.arguments = np.array([bstring(input_node1), bstring(input_node2)])
    dgrp._v_attrs.dim1 = dim1
    create_array(h5, dgrp, 'id', obj=ids)
    return name

def AngleCoord(h5, nodename, input_node, ids):
    potential = h5.root.input.potential
    base = 'Angle'
    if nodename is None:
        name = base
    else:
        name = '{}_{}'.format(base, nodename)
    agrp = h5.create_group(potential, name)
    agrp._v_attrs.arguments = np.array([bstring(input_node)])
    create_array(h5, agrp, 'id', obj=ids)
    return name

def TorsionCoord(h5, nodename, input_node, ids):
    potential = h5.root.input.potential
    base = 'Dihedral'
    if nodename is None:
        name = base
    else:
        name = '{}_{}'.format(base, nodename)
    dgrp = h5.create_group(potential, name)
    dgrp._v_attrs.arguments = np.array([bstring(input_node)])
    create_array(h5, dgrp, 'id', obj=ids)
    return name

def COMCoord(h5, nodename, input_node, n_dim, index_pos, border, index_dim):
    potential = h5.root.input.potential
    base = 'GroupCenter'
    if nodename is None:
        name = base
    else:
        name = '{}_{}'.format(base, nodename)
    grp = h5.create_group(potential, name)
    grp._v_attrs.arguments = np.array([bstring(input_node)])
    grp._v_attrs.n_dim     = n_dim
    grp._v_attrs.n_group   = len(border) - 1
    assert index_pos.size == border[-1]
    create_array(h5, grp, 'index_pos', obj=index_pos)
    create_array(h5, grp, 'border',    obj=border)
    create_array(h5, grp, 'index_dim', obj=index_dim)
    return name

def Const3DCoord1(h5, nodename, input_node, values):
    potential = h5.root.input.potential
    base = 'Const3D'
    if nodename is None:
        name = base
    else:
        name = '{}_{}'.format(base, nodename)
    grp = h5.create_group(potential, name)
    grp._v_attrs.arguments = np.array([bstring(input_node)])

    values = np.array(values)
    assert values.ndim == 2
    ids = np.arange(values.shape[0])
    create_array(h5, grp, 'id', obj=ids)
    create_array(h5, grp, 'value', obj=values)

    return name

def MovingConst3DCoord1(h5, nodename, input_node, time_initial, time_step, velocities, values):

    potential = h5.root.input.potential
    base = 'MovingConst3D'
    if nodename is None:
        name = base
    else:
        name = '{}_{}'.format(base, nodename)

    grp                               = h5.create_group(potential, name)
    grp._v_attrs.arguments            = np.array([bstring(input_node)])
    grp._v_attrs.time_initial         = time_initial
    grp._v_attrs.time_step            = time_step
    grp._v_attrs.step_initial         = 0   # time steps before this run; a continuation resumes /output/tip_step

    ids = np.arange(velocities.shape[0])
    create_array(h5, grp, 'id', obj=ids)
    create_array(h5, grp, 'velocities', obj=velocities)
    create_array(h5, grp, 'start_pos', obj=values)

    return name

def SpringPotential(h5, nodename, input_node, dim, is_pbc, box_len, ids, equil_dist, spring_const):
    potential = h5.root.input.potential
    base = 'Spring'
    if nodename is None:
        name = base
    else:
        name = '{}_{}'.format(base, nodename)
    grp = h5.create_group(potential, name)
    grp._v_attrs.arguments = np.array([bstring(input_node)])
    grp._v_attrs.dim1      = dim
    grp._v_attrs.pbc       = is_pbc
    grp._v_attrs.box_len   = box_len
    create_array(h5, grp, 'id',           obj=ids)
    create_array(h5, grp, 'equil_dist',   obj=equil_dist)
    create_array(h5, grp, 'spring_const', obj=spring_const)
    return name

def WallSpringPotential(h5, nodename, input_node, dim, ids, equil_dist, spring_const, wall_type):
    potential = h5.root.input.potential
    base = 'WallSpring'
    if nodename is None:
        name = base
    else:
        name = '{}_{}'.format(base, nodename)
    grp = h5.create_group(potential, name)
    grp._v_attrs.arguments = np.array([bstring(input_node)])
    grp._v_attrs.dim1      = dim
    create_array(h5, grp, 'id',           obj=ids)
    create_array(h5, grp, 'equil_dist',   obj=equil_dist)
    create_array(h5, grp, 'spring_const', obj=spring_const)
    create_array(h5, grp, 'wall_type',    obj=wall_type)
    return name
