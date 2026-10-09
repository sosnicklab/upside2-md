# Progress from the MacBook Pro session of 2026-10-09 (to be merged into progress.md)

Separate from `progress.md` while the Mac Studio owns it (user, 10-09). No tracked file was edited;
`git status` is clean.

* **11:52, remote status, read only.** Mac Studio owns the jobs (last pass 11:47). Reported squeue and
  WATCH_STATUS: ff30_gdepth_dt009 converged (local ff_3.0); bio runs at steps 57-64 of 76; 32
  ff_3.0_gdepth benchmark arms running, 4 glpG chains queued.
* **11:55-12:30, glycine map figures.** `~/Downloads/gly_rama_ff30.png` (format of the earlier two),
  then all three maps on one zero (`*_nonhelix_zero.png`) and a US Letter print page
  `~/Downloads/gly_rama_nonhelix_zero_page.pdf` (`scratchpad/plot_gly_rama_ff30.py`,
  `plot_gly_rama_nonhelix_zero.py`). findings_macbookpro M1.
* **11:55-13:05, local glpG pair.** ff2.1 and ff3.0 patched (`patch_glpg.py` md5 70589119...) and run,
  seed 1, TM4-test flags; bitwise equal to the cluster's seed-1 runs. Figures
  `~/Downloads/glpg_tm4_ss_ff21_vs_ff30.png` (slide-13 format, phi/psi) and
  `glpg_tm4_dssp_ff21_vs_ff30.png` (DSSP alpha), six TM helices marked from our PDB. M2, M3.
* **13:10-13:40, TM helices across sets.** Backbones of 9 sets x 12 seeds streamed from `runs_cov`
  (no compute on the login node beyond reading, ulimit 4 GB): `sets_bb*.npz`, `sets_full.npz`. TM
  breakdown, attribution and H-bond tables. M4.
* **13:40-14:10, TM1 diagnosis.** Lipid environment and energy terms at the opening, time-matched over
  12 seeds; Upside's membrane H-bond term evaluated on stored frames. M5. User chose option A.
* **13:45-13:52, option A built and verified** (`scratchpad/glpg_tm1_hbmem/`): energy identity and
  finite-difference forces pass. 24 runs started 13:48; the ff3.0 twelve stopped at t 60 for the
  16:00 shutdown (user), ff2.1 twelve continue (~15:35).
* **13:50-14:05, transient backbone excursions found** in the baseline sets (M6); dense-frame replay
  of ff2.1 seed 1 with pair diagnostics started.
* **14:20-14:35, VTF and the phantom O.** `~/Downloads/glpG_79HIS_ff3.0_T080_s1.vtf` (ff3.0, seed 1, no
  TM1 term). The stationary atom the user saw is residue 210's unused O slot (M7); fixed in
  `py/martini_extract_vtf.py` (since committed by the user, 14a27cc3) and the VTF regenerated.
* **15:39-15:50, option A result.** ff2.1 + term finished (12 seeds, 401 frames each). Helicity up in
  every TM helix (TM1 0.79 -> 0.90, TM4 primary 0.609 -> 0.788), but 5 seeds jump > 3000 E_up and two
  tear TM4's backbone (C-N 9.5 A). M8. `analyze_hbmem.py` limited to the finished set and its frame
  check tightened (401 frames, t 0); output `analyze_hbmem.txt`. Positions saved compactly as
  `runs_hbmem_ff21_pos.npz` and copied to the cluster with the logs, kick outputs and these records.
* **15:40, DSSP figure question.** The DSSP-alpha figure differs from slide 13 at t 0 because of the
  criterion (31 residues, mostly 3-10 and turns outside the TMs, plus TM2's pi label); the phi/psi
  figure matches the slide's first frame on 208/208 residues.
