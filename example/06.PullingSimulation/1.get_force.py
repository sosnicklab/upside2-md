import sys,os
import numpy as np
import tables as tb
import matplotlib.pyplot as plt

def _output_groups(t):
    i=0
    while 'output_previous_%i'%i in t.root:
        yield t.get_node('/output_previous_%i'%i)
        i += 1
    if 'output' in t.root:
        yield t.get_node('/output')
        i += 1

upside_to_pN   = 41.4

#-------------------------------------

AFM_table  = sys.argv[1]

pdb_id     = '1qhj'
sim_id     = 'pulling_test'
n_rep      = 8 
work_dir   = './'

#-------------------------------------
output_dir = "{}/outputs".format(work_dir)
result_dir = "{}/results".format(work_dir)
run_dir    = "{}/{}".format(output_dir, sim_id) 

if not os.path.exists(result_dir):
    os.makedirs(result_dir)

h5_files   = []
for j in range(n_rep): 
    h5_file = "{}/{}.run.{}.up".format(run_dir, pdb_id, j)
    h5_files.append(h5_file)

#-------------------------------------

fields = [ln.split() for ln in open(AFM_table)]
header1 = 'residue spring_const tip_pos_x tip_pos_y tip_pos_z pulling_vel_x pulling_vel_y pulling_vel_z'
header2 = 'residue spring_const pulling_vel_x pulling_vel_y pulling_vel_z'
actual_header = [x.lower() for x in fields[0]]

fields = fields[1:]
n_spring = len(fields)

atom = []
k    = []
v    = []
for i,f in enumerate(fields):
    atom.append( int(f[0])*3 + 1 )
    k.append( float(f[1]) )
    if len(fields[0]) == 8:
        v.append( [float(x) for x in (f[5],f[6],f[7])] )
    else:
        v.append( [float(x) for x in (f[2],f[3],f[4])] )
atom = np.array(atom)
v    = np.array(v)
k    = np.array(k)

ndx = np.where(v!=0)

dim_name = ['x', 'y', 'z']

for ii in range(ndx[0].size):
    aid = ndx[0][ii]
    aa  = atom[aid]
    dim = ndx[1][ii]

    for rid, h5 in enumerate(h5_files):
        with tb.open_file(h5, 'r') as t:
            for g_no,g in enumerate(_output_groups(t)):
                if g_no == 0:
                    pos     = g.pos[:,0, aa, dim]
                    tip_pos = g.tip_pos[:, aid, dim]
                else:
                    pos     = np.concatenate([ pos,     g.pos[1:,0,aa, dim]   ])
                    tip_pos = np.concatenate([ tip_pos, g.tip_pos[1:, aid, dim] ])

        force = k[aid] * (tip_pos-pos) * upside_to_pN
        start = pos[0]*1
        extension = pos - start
        tip_exten = tip_pos - start
        fname = '{}/{}-{}-replica{}-res{}-{}.dat'.format(result_dir, pdb_id, sim_id, rid, (aa-1)//3, dim_name[dim] )
        with open(fname, 'w') as out:
            out.write('# Extension (nm) Force (pN) TipPos (nm)\n')

            for i,e in enumerate(extension):
                out.write('{} {} {}\n'.format(e, force[i], tip_exten[i]))

