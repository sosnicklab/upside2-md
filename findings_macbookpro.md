# Findings from the MacBook Pro session of 2026-10-09 (to be merged into findings.md)

Written on the MacBook Pro while the Mac Studio owns `findings.md` (user, 10-09: keep this session's
records in separate files; another session merges them). Scripts, outputs and data are in
`scratchpad/glpg_ff21_vs_ff30/` unless named otherwise.

## M1. ff3.0's glycine map is NDRD plus one depth pair (checked against the file)

- `parameters/ff_3.0/rama.dat` (md5-identical to `$P/parameters/ff_3.0_gdepth`) differs from
  `common/rama.dat` only in the 37 GLY|X coil maps; the sheet group, every other map and the pooled
  GLY|ALL entry are NDRD's.
- All 37 carry one identical shift (to 1.2e-6). Fitted to the trainer's own basin weights
  (`ff30_gdepth_dt009/trainer/rama_basin.py`) it is alpha_R -0.1068, alpha_L +0.5169, the round-4
  offsets in `rama_rounds.txt`; the trainer's write rule reproduces every stored map to 1e-6.
- **A label value is not a depth.** The library renormalises each map and the config writer subtracts
  the Boltzmann-weighted mean (`upside_config.py:841`), so raising alpha_L moves every other cell by a
  constant: -0.155 from renormalisation and -0.167 from the mean, -0.322 in all. On the 2,492 cells the
  offsets do not touch, ff3.0 minus ff2.1 is -0.321 to -0.322. PII - alpha_L moves by exactly -0.517
  and PII - alpha_R by +0.107.
- For figures, zero each map at the Boltzmann-weighted mean outside the trainer's helical basins
  (helix weight < 1e-4, 2,072 cells): ff3.0 then equals ff2.1 there to 3e-5. On that zero, alpha_L /
  alpha_R / PII are -3.0 / -1.6 / -0.6 (ff2.1), -2.5 / -1.7 / -0.6 (ff3.0), -1.2 / -0.9 / -0.9
  (BioEmu-fitted). `scratchpad/plot_gly_rama_nonhelix_zero.py`.

## M2. The local glpG pair reproduces the TM4-test sets frame for frame

ff2.1 and ff3.0 patched into the live 79HIS seed with the fixed `patch_glpg.py`, run with
`run_glpg.sh`'s flags (T 0.80, 4000 tu, seed 1) on the MacBook Pro: all 401 log frames are identical
to `runs_cov/79HIS_ff21_released_T080_s1.log` and `79HIS_d9_03_T080_s1.log`. KE/1.5kT 1.007 and 1.002;
no total-potential jump above 3000. Seed 1 is ff3.0's 10th-best of 12 by last-block DSSP alpha and
ff2.1's 6th, so the pair shows a larger contrast than the 12-seed sets (0.726 against 0.610, p 0.25).

## M3. glpG's six TM helices, from the structure we simulate

`scratchpad/local_popg_79HIS/glpG-RKRK-79HIS.pdb` (= `pdb_staging/glpG-RKRK-79HIS.pdb`, md5
`e8c3e1a5...`) is the live seed's t = 0 protein (CA RMSD 0.00 A). Its DSSP helices: TM1 29-48, TM2
82-103 (DSSP I at 92-96), TM3 105-127, TM4 135-151, TM5 161-175, TM6 185-207, as findings 11 gives.
CA depth spans 21-30 A against PO4 planes at +-21 A. 19-26 and 50-57 lie flat in the interfaces
(4 and 7 A spans) and are not TM. Its B-factor column runs 88-96 in TM2 and 39-58 at the termini,
the reverse of crystallographic B-factors; it reads like a prediction confidence (unconfirmed).

## M4. Which TM helices lose structure, and from what (12 seeds per set, last block, DSSP)

`tm_attribution.txt`, `tm1_tm2_hbonds.txt`; any-helix = DSSP H, G or I.
- **TM1 loses structure in every set** (any-helix 0.75-0.85; ff2.1 0.79, ff3.0 0.85, control c9_01
  0.75) and is still falling at the end (ff2.1 0.89, 0.85, 0.80, 0.80 by block; ff3.0 0.94, 0.92,
  0.87, 0.85, alpha 0.83 -> 0.69). Its middle, 32-38, sits at the bilayer centre (CA z -7 to +3 A).
  ff2.1 unwinds it (about 39% of last-block i->i+4 H-bonds broken, 5% shifted to i->i+5); ff3.0
  turns more of it into a pi-bulge (20% i->i+5, about 19% broken). TM1 has no glycine.
- **TM2 is lower in ff3.0 than ff2.1** (any-helix 0.79 against 0.91, p 0.04; alpha p 0.08), at its
  N-terminal end 82-85 (CA z +3 to +6 A) and at 91-96. The 82-85 loss comes with training of the
  non-glycine terms: gdepth_start (ff2.1 + the glycine map) is 0.86, the control c9_01 (ff2.1's own
  glycine map, trained otherwise) 0.79 (p 0.01), and the frozen-H-bond twins lean higher (bz_02 0.88
  against b9_02 0.80; dz_00 against d9_00 +0.13 at 82-85). None of the twin differences resolves.
- **TM2's start has no gap.** At 91-96 every i->i+4 O...N is 2.8-3.5 A and phi/psi are helical; DSSP
  calls 92-96 pi from one bifurcated O92...N97 contact at 3.44 A. ff2.1 keeps the i->i+4 bonds
  (0.72 of acceptor-frames); ff3.0 loses them (0.54, about 29% broken; any-helix 0.95 -> 0.75 by
  block, 6 of 12 seeds below 0.80). Only d9_03 does this: c9_01, b9_02, bz_02, gdepth_start keep
  91-96 at 0.89-0.97, and d9_02 was at 0.87. Not resolved against ff2.1 (p 0.16).
- **TM5** is within noise of ff2.1 (0.88 against 0.92, p 0.34).
- 12-seed TM4 values reproduce the cluster tables (d9_03 0.726, ff21_released 0.609-0.610,
  gdepth_start 0.657), which checks the reader.

## M5. TM1's middle opens because nothing charges for an unpaired backbone in the acyl core

`tm1_diag.txt`, `tm1_sets.txt`, `tm1_hbmem.txt`.
- **Lipid gain is not the driver.** Time-matched (last block, 12 seeds), TM1's middle has the same
  environment intact or open: ff2.1 about 18 tail beads within 7 A either way and BB-env energy -43.7
  against -44.3 E_up (hb4 against BB-env r +0.06); ff3.0 -40.6 against -44.0. No headgroup or ion
  comes near. The single-seed contrast (-25.8 against -41.6) was tails accumulating with time.
- **What opening costs in the hybrid:** only Upside's H-bond energy, E_alpha -1.96 E_up per bond
  (about +8.4 E_up for the 4.3 bonds the middle loses). MARTINI's BB typing is fixed at the starting
  secondary structure (TM1: 17 N0, 3 C5), so an opened backbone keeps its helix-typed lipid
  attraction (N0-C1 well -1.18 E_up; MARTINI's coil type P5 has -0.17).
- **Upside's own membrane model charges for exactly this.** `hb_membrane_potential`
  (`src/membrane_potential.cpp:488`, `parameters/ff_2.1/membrane.h5` hb_energy) gives each backbone
  donor and acceptor (1-p)^2 f_unpaired(z) + (1-(1-p)^2) f_paired(z). f_unpaired - f_paired is about
  +2.2 (donor) and +1.8 (acceptor) E_up within ~5 A of the midplane and slightly negative in water.
  Its state-dependent part, evaluated on the stored frames with protein_hbond's p (0.74 intact, 0.35
  open in ff2.1), would add +7.0 / +8.4 / +10.3 E_up to opening 32-38 at half-thickness 12.7 / 14.0 /
  15.9 A (ff3.0 +4.4 / +5.3 / +6.4). The hybrid has no such term: it charges about half of what
  Upside's membrane model says opening TM1's middle costs.
- The bilayer centre sits within 1 A of z = 0 and drifts 1-2 A over 4000 tu (PO4 mid-plane).
- Tested with ff2.1 (M8): the term raises TM1's helicity but tears TM4's backbone in 2 of 12 seeds.

## M6. Transient backbone excursions sit where the helices fail (2026-10-09, 13:50-14:05)

`scratchpad/glpg_tm1_hbmem/cn_events_vs_tm1.txt`; ff21_released and d9_03, 24 seeds, all frames.
- 61 of 9,624 stored frames hold a peptide C-N above 2 A (up to 6.0 A; median C-N 1.327, normal
  frames' maximum 1.6-1.8). Each lasts one stored frame (frames are 370 steps apart) and recovers.
  They cluster in TM1 33-36, TM4 137-146 and the 74-82 loop, several consecutive bonds at once.
- ff2.1 seed 1, frame 41 (t 410): C-N 5.20 A at 34; Spring_bond 277 -> 1105, Spring_angle 259 -> 738,
  backbone_pairs 6 -> 50, protein potential +1477 E_up, all back at frame 42; martini_potential not
  unusual. This is the frame TM1's middle first opens. The jump is under the 3000 threshold of the
  TM4 test's jump scan, so these events are not in its counts.
- 16 excursions in TM1's middle: its i->i+4 H-bonds average 5.10 over the 5 frames before and 3.83
  over the 5 after, mostly 1-3 of 6 at the event frame. Of 7 sustained openings, 2 have an excursion
  within the window, a lower bound at this frame spacing.
- **Not a MARTINI core collision at the snapshot:** closest BB-environment approach at the excursion
  frames is median 5.09 A (min 4.10) against 5.15 A three frames earlier (healthy 4.03, blow-up 2.43,
  findings 2.4). The strain is in Upside's own bonded terms.
- **The dense-frame replay localises it** (`scratchpad/glpg_tm1_hbmem/kick/`: `replay.sh`,
  `analyze_replay.txt`, `force_by_node.txt`). ff2.1 seed 1 replayed to t 420 with frames every 10
  steps reproduces the stored run bit for bit (43 shared frames, max difference 0). The onset is
  between t 407.70 and 407.97: MET34's backbone moves 3.6 A in 10 steps against 0.5-0.9 A for every
  other atom, then up to 9.5 A per 10 steps for 0.5 tu, C-N 7-9 A, Spring_bond up to 9,067 E_up
  (normal 220-290); it settles by t 415 with TM1's middle open.
- **No MARTINI pair spike precedes it.** `UPSIDE_MARTINI_PAIR_DIAG` (thresholds 3.6 A, 200 E_up/A)
  prints nothing between t 395 and 409; its one report, 3.43 A and 676 E_up/A at t ~409, comes after
  the backbone has flown apart, and is far from the blow-up regime (1.3e5 E_up/A at 2.85 A).
- **No term shows a precursor at the stored frames.** With each leaf potential node removed from a
  copy of the input, every term's force on N/CA/C of 33-35 is at most ~47 E_up/A through t 407.70;
  by 407.97 Spring_bond is restoring the torn geometry (137, then 437 E_up/A). The only unusual
  values are two backbone_pairs contacts on N33 (47 and 37 E_up/A at t 406.08 and 407.16), zero
  otherwise. The impulse lies inside those 10 steps (~40 force evaluations) and is not seen at
  either frame. Not tested: the side-chain/lipid 1-body table separately (it sits inside rotamer),
  and the steps themselves (a replay saving every step takes ~2.6 h here at ~0.6 s per frame).

## M7. The VTF writer drew a phantom atom: the C-terminal O slot (fixed 2026-10-09)

The O slot of glpG residue 210 never moves (0.0000 A over 401 frames, 19 A from its C at the end). The
engine writes a residue's O slot only where infer_H_O has an acceptor on that residue's carbonyl C
(`src/martini_hybrid.cpp`, `resolve_bb_o_hbond_elements`); the C-terminal residue has none, so the slot
keeps its input coordinate. Nothing reads it: it is in none of the 8.4 M MARTINI pairs (only BB beads
are), no potential node references it, and it is not in `/input/brownian` (no O slot is). It has no
effect on the dynamics or on any analysis here (TM6 ends at 207). `py/martini_extract_vtf.py` wrote it
as an atom. Fixed there by the engine's own rule: `build_backbone_projection_map` records
`carrier_present` (N, CA, C always; O where an acceptor sits on the residue's C), and the atom list
(`backbone_output_atoms`, shared by modes 1 and 2), the frame assembly and the bonds (`backbone_bonds`,
which replaces `mode2_backbone_bonds`) follow it. glpG now has 839 backbone atoms (209 O), bonds
N-CA 210, CA-C 210, C-O 209, C-N 209; O positions equal the trajectory's to 0.001 A; modes 1 and 2
both run. In the user's commit 14a27cc3.

## M8. The membrane H-bond term holds TM1 but tears TM4 (ff2.1, 12 seeds, 2026-10-09 15:40)

`scratchpad/glpg_tm1_hbmem/analyze_hbmem.txt`; ff2.1 + term (option A, half-thickness 15.9 A) against
ff21_released, last block, DSSP any-helix (alpha), two-sided Mann-Whitney on 12 seeds each.
- **Every TM helix is more helical:** TM1 0.79 -> 0.90 (p 0.03), TM4 0.71 -> 0.82 (p 0.02), TM5
  0.92 -> 0.97 (p < 0.01), TM2/TM3/TM6 +0.01 to +0.04 (n.s.); TM4 primary 0.609 -> 0.788 (p 0.03).
  Seeds with TM1 below 0.80 fall from 6 to 2; TM1 32-38 i->i+4 H-bonds 0.56 -> 0.66.
- **It fails the health check.** Protein potential jumps by 3,800-17,900 E_up between stored frames in
  5 of 12 seeds (s1, s2, s9, s11, s12; the local no-term ff2.1 seed 1 has none above 2,500). In s1
  (t 1680) and s2 (t 2940) TM4's backbone tears over 137-148, C-N up to 9.5 A, for 2-3 stored frames;
  the others are single-bond tears at 36, 198 and 74-84. Frames with any C-N > 2 A: 27 against 28,
  but 8 above 4 A against 4, and no TM4 tear of that size in ff21_released. Rg and the protein's
  depth are unchanged; KE/1.5kT 1.001-1.073 (s2 high, from its tear).
- Cause of the tears not measured. The term acts through protein_hbond's p on N, CA and C and its
  forces matched finite differences on stored frames (verify_hbmem.txt), so a wrong derivative is
  not the explanation; whether it raises the stored strain that M6's excursions release is untested.
- Not run: ff3.0 + term (stopped at t 60 for the shutdown).

## Lessons (this session)

- **An energy-structure correlation along one trajectory is confounded with time.** In ff2.1 seed 1,
  "opening is rewarded by lipid" came from tails accumulating while the run aged; time-matched
  frames over 12 seeds showed no difference. Compare states at matched times across seeds before
  naming a driver.
- **Annotate from the structure we simulate, never from online data** (user, 10-09: our construct has
  no online counterpart). Derive helix ranges from DSSP on the seed's own PDB and check them in the
  script.
- **When another computer owns the shared md files, write this session's records to separate files**
  (user, 10-09); git must show `plan.md`, `findings.md`, `progress.md`, `remote_jobs.md` unchanged.
- **A DSSP pi label is not a gap.** Check the i->i+4 and i->i+5 O...N distances before reading a
  white band in a DSSP-alpha figure as a helix break.
- **Compare a figure with a reference only under the reference's criterion.** The slide used a phi/psi
  box; a DSSP-alpha companion differed from frame 0 on 31 residues (3-10, turns, the TM2 pi label)
  and the user read it as a different structure. Name the criterion in the companion's title.
- **An interim health count is not a trend.** Mid-run, the term set had fewer C-N excursions; the
  finished set has the largest tears of any. Report health only from finished runs.
- SciencePlots' `science` style sets `savefig.bbox: tight`, which crops a fixed page size; set
  `'savefig.bbox': 'standard'` for a print page.
