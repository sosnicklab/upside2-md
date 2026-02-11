# Findings

No external findings logged yet.
# Gromacs CHARMM-GUI MARTINI (1ubq) .mdp settings (Downloads/charmm-gui-7077081595/gromacs)

## Minimization
- step4.0_minimization.mdp: integrator=steep, dt=0.020, nsteps=3000, tcoupl=v-rescale (tau_t=1.0, ref_t=303.15), Pcoupl=berendsen (tau_p=5.0, compressibility=4.5e-5, ref_p=1.0), reaction-field (epsilon_r=15, rcoulomb=1.1, rvdw=1.1), soft-core enabled (free-energy yes, init-lambda=0.01, sc-alpha=4, sc-power=2, sc-coul yes).
- step4.1_minimization.mdp: integrator=steep, dt=0.020, nsteps=3000, same TC/P settings, no NORMANG define, no soft-core block.

## Equilibration/Production
- step4.2_equilibration.mdp: integrator=md, dt=0.020, nsteps=50000, Pcoupl=berendsen (tau_p=5.0, compressibility=4.5e-5).
- step5_production.mdp: integrator=md, dt=0.020, nsteps=1000000, Pcoupl=Parrinello-Rahman (tau_p=12.0, compressibility=4.5e-5).

## Unit conversion reference (Upside)
- 4.5e-5 bar^-1 compressibility corresponds to 2.1782 Ang^3/E_up (from 3e-4 bar^-1 = 14.521180763676 Ang^3/E_up scaling).
