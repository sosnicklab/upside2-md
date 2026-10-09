# Plan from the MacBook Pro session of 2026-10-09 (to be merged into plan.md)

Separate from `plan.md` while the Mac Studio owns the shared records (user, 10-09).

## Goal

Hold glpG's TM1 in the hybrid without retraining the force field (user, 10-09), by correct physics:
no scaling or disabling of BB-env or SC-env, no tuned parameters, no guards.

## Diagnosis (done; findings_macbookpro M4-M5)

TM1's middle (32-38) sits at the bilayer centre and loses structure in every force field, ff2.1 included.
Lipid gain is not the driver. Opening it costs only Upside's soluble-calibrated H-bond energy: the
hybrid has no cost for an unpaired backbone NH or CO in the acyl core, and its MARTINI backbone typing
is fixed at the starting secondary structure. Upside's own implicit-membrane model has that cost
(`hb_membrane_potential`), and on the stored frames it would add +7 to +10 E_up to the opening.

## Approach (user chose A, 2026-10-09 13:45)

**A. Add Upside's membrane backbone H-bond term, its H-bond-state-dependent part only.**
- E = sum over backbone donors and acceptors of (1-p)^2 [f_unpaired(z) - f_paired(z)], p from the
  hybrid's existing `protein_hbond`, f from `parameters/ff_2.1/membrane.h5` hb_energy unchanged. The
  paired baseline is left out because MARTINI's helix-typed BB already describes a paired backbone
  against lipid.
- First test is config-only: the engine's existing `hb_membrane_potential` node, written into a copy
  of the patched glpG input with coeff = [hb_energy[t,0] - hb_energy[t,1], 0], use_curvature 0, the
  bilayer centred at z = 0 (it stays within 1-2 A over 4000 tu). No C++ change, master parity untouched.
- Half-thickness from our own bilayer, not borrowed: the glycerol-ester (GL1/GL2) mid-plane, the
  hydrocarbon-core boundary Upside's membrane thickness denotes; the C1 plane (12.4 A) would give a
  smaller penalty (tm1_hbmem.txt has the sensitivity).
- Production would need the reference to follow the bilayer's centre (a node reading the lipid
  mid-plane) and a builder option in `py/martini_prepare_system*.py`; decided after the test.

**B. Make each BB bead's MARTINI type follow its H-bond state** (N0 paired, P5 unpaired). MARTINI's
own endpoints, but a new C++ coupling in `martini_potential` with derivatives through `protein_hbond`;
on the stored frames it would charge about +38 E_up for the opening, far more than A.

## Execution phases (A)

- [x] Half-thickness 15.9 A: mid-plane of GL1 (16.42 A) and GL2 (15.40 A) over the second half of the
      local ff2.1 run; the first acyl beads sit at 12.1-12.7 A (seed: GL1 17.6, GL2 16.3).
- [x] Term built into copies of both inputs (`scratchpad/glpg_tm1_hbmem/add_hbmem.py`). Verified on
      stored frames (`verify_hbmem.txt`): E_new - E_old equals the node value and an independent
      evaluation to 1e-3 E_up; node forces match finite differences to ~1e-3. +86.5 E_up at t 0.
- [x] Test, ff2.1 + term, seeds 1-12 (finished 15:39; findings_macbookpro M8): every TM helix more
      helical (TM1 +0.11, TM4 primary +0.18), but TM4's backbone tears in 2 seeds and the protein
      potential jumps > 3000 E_up in 5. Not acceptable as it stands.
- [ ] Find why the term produces the tears (replay s1 around t 1680 or s2 around t 2940 with dense
      frames and the per-node force ablation of kick/force_by_node.py) before any further use of it.
- [ ] ff3.0 + term, 12 seeds, on another computer or the cluster, after the tears are explained.
- [ ] Report; decide production implementation with the user.

## Second lead: transient backbone excursions (findings_macbookpro M6)

- [x] Dense-frame replay of ff2.1 seed 1 to t 420 (`kick/replay.sh`): bitwise reproduction; onset
      between t 407.70 and 407.97 at MET34's backbone; no MARTINI pair spike before it; no term's force
      on 33-35 abnormal at the stored frames (`kick/force_by_node.txt`). findings_macbookpro M6.
- [ ] Replay saving every step from t ~407.6 (whole run from t 0, ~2.6 h at the current output cost,
      or with a cheaper output path), then the per-term force decomposition of each step, with the
      side-chain/lipid 1-body table split out of rotamer. Next computer.

## Known errors / blockers

- The membrane model's f was trained in Upside's implicit membrane, with its own burial terms; using
  its state-dependent part in the hybrid is a modelling decision for the user (physics change).
- This computer is in stand-by for the cluster jobs (remote_jobs.md Handoff); local runs here do not
  touch the TM4 queue.
- The MacBook Pro shuts down ~16:00 (user); anything unfinished then moves to another computer.
- The added term is anchored at z = 0 and acts on the protein only; watch the protein's depth against
  the explicit bilayer in the test (`analyze_hbmem.py` health block).
