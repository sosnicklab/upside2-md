import sys, os
import numpy as np
import tables as tb
import mdtraj as md

upside_path = os.environ['UPSIDE_HOME']
upside_utils_dir = os.path.expanduser(upside_path+"/py")
sys.path.insert(0, upside_utils_dir)

import mdtraj_upside as mu
import upside_engine as ue
import matplotlib.pyplot as plt

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('top_h5',            help='Input top file (.up or .h5 file)')
    parser.add_argument('input_h5',          help='Input simulation file')
    parser.add_argument('output_npy',        help='Output npy file')
    parser.add_argument('--use-TM-region',   default=False, action='store_true', help='(default false) treat the NH on the protein-lipid surface as protected')
    parser.add_argument('--report-raw-data', default=False, action='store_true', help='(default false) report the raw data (H-bond score, burial level, prot-lipid surface score)')
    parser.add_argument('--stride',          type=int,   default=1, help='(default 1) Stride for reading file')
    parser.add_argument('--start',           type=int,   default=0, help='(default 0) Initial frame')
    parser.add_argument('--residue',         type=str,   default=None, help='(default none) the file used to store the residue id')
    parser.add_argument('--criterion1',      type=float, default=0.01, help='(default 0.01) to judge whether NH is H-bonded. bigger than the criterion means H-bonded')
    parser.add_argument('--criterion2',      type=float, default=0.05, help='(default 0.05) to judge whether NH is H-bonded by side chain aceptor. bigger than the criterion means H-bonded')
    parser.add_argument('--criterion3',      type=float, default=5.00, help='(default 5.00) to judge whether NH is exposed. bigger than the criterion means buried')
    parser.add_argument('--criterion4',      type=float, default=0.50, help='(default 0.50) to judge whether NH is exposed to the lipid. bigger than the criterion means exposed to the lipid')
    args = parser.parse_args()

    engine = ue.Upside(args.top_h5)

    with tb.open_file(args.top_h5, 'r') as t:
        donor = t.root.input.potential.infer_H_O.donors.residue[:]
        n_donor = donor.size
        # for side chain buiral level
        weight = t.root.input.potential.sigmoid_coupling_environment.weights[:20]
        # for side chain-NH H-bond
        weight2 = weight*0
        weight2[3] = 1. # ASP
        weight2[6] = 1. # GLU

    traj    = mu.load_upside_traj(args.input_h5, top=args.top_h5)
    bb      = traj.top.select('name C or name N or name CA')
    traj_bb = traj.atom_slice(bb)[args.start::args.stride]
    N = traj_bb.n_frames

    print ("{} frames are used".format(N))

    Hbond1 = []
    Hbond2 = []
    Burial = []
    Surf   = []
    for i in range(N):
        p = engine.energy(traj_bb.xyz[i]*10)

        hb  = engine.get_output('protein_hbond')[:,6]
        Hbond1.append(hb[:n_donor])
    
        # BL from backbone atoms
        bl1 = engine.get_output('hbbb_coverage')[:,0] 

        # BL from side chain beads
        burial_level2 = engine.get_output('environment_coverage_hb')[:,0] 
        n_bb = burial_level2.size/20
        bl2 = np.zeros(n_donor)
        for ib in range(n_donor):
            for aa in range(20):
                bl2[ib] += weight[aa] * burial_level2[ib*20+aa]
        Burial.append(bl1+4.6*bl2)

        # BL from  ASP or GLU (to evaluate H-bonds from side chain acceptors)
        bl3 = np.zeros(n_donor)
        for ib in range(n_donor):
            for aa in range(20):
                bl3[ib] += weight2[aa] * burial_level2[ib*20+aa]
        Hbond2.append(bl3)

        if args.use_TM_region:
            surf = engine.get_output("surface")[donor,0]
            Surf.append(surf)

    Hbond1 = np.array(Hbond1)
    Hbond2 = np.array(Hbond2)
    Burial = np.array(Burial)

    HB1 = Hbond1*0.
    HB1[Hbond1>args.criterion1] = 1.
    HB2 = Hbond2*0.
    HB2[Hbond2>args.criterion2] = 1.
    BL = Burial*0.
    BL[Burial>args.criterion3] = 1.

    PS = HB1 + HB2 + BL

    if args.use_TM_region:
        Surf = np.array(Surf)
        Su = Su*0.
        Su[Surf>args.criterion4] = 1.
        PS += Su

    PS[PS>1.] = 1.

    if '.npy' in args.output_npy:
        output_base = args.output_npy[:-4]
    else:
        output_base = args.output_npy[:]

    np.save(output_base, PS)

    if args.report_raw_data:
        np.save('{}_Hbond.npy'.format(output_base), Hbond1)
        np.save('{}_CTerHbond.npy'.format(output_base), Hbond2)
        np.save('{}_BurialLevel.npy'.format(output_base), Burial)
        if args.use_TM_region:
            np.save('{}_ProtLipidSurf.npy'.format(output_base), Surf)

    if args.residue:
        np.savetxt(args.residue, donor, fmt='%i')

if __name__ == '__main__':
    main()
