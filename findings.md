# Findings

- 2026-02-12: Dry MARTINI Gromacs workflow files are present at `/Users/yinhan/Downloads/charmm-gui-7090685331/gromacs` with sequence: `step6.0_minimization.mdp`, `step6.1_minimization.mdp`, `step6.2` to `step6.6` equilibration, and `step7_production.mdp`.
- 2026-02-12: Existing bilayer workflow script already performs staged minimization/equilibration/production but was tuned for prior MARTINI assumptions.
- 2026-02-12: Dry MARTINI mdp settings indicate semi-isotropic coupling for relaxation (`tau-p=4.0`, membrane-plane compressibility equivalent to `3e-4 bar^-1`, normal-axis fixed for implicit systems) and a pressure-coupling strategy that transitions from Berendsen to Parrinello-Rahman as described in the paper note.
- 2026-02-12: UPSIDE barostat previously used a single compressibility scalar for all axes; true dry MARTINI membrane coupling requires axis-specific compressibility (`xy` vs `z`).
- 2026-02-12: `ff_dry` uses dry MARTINI v2.1 file names (`dry_martini_v2.1.itp`, `dry_martini_v2.1_lipids.itp`, `dry_martini_v2.1_ions.itp`) and macro-based aliases in `[ bonds ]`/`[ angles ]`, so parser support for `#define`-backed parameters is required.
- 2026-02-12: Dry lipid ITP (`dry_martini_v2.1_lipids.itp`) relies on preprocessor branches (`#ifndef EXP_DOPC`/`#else`), and only one branch should be parsed.
