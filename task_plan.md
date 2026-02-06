# Plan: MARTINI Water Diffusion Scan on Slurm

Project Goal: Update the MARTINI water diffusion scan scripts and workflow to run on Slurm, sweep temperature and THERMOSTAT_TIMESCALE ranges, and generate diffusion plots based on the last 50% of production.

## Architecture & Key Decisions
- Keep the scan workflow in `example/16.MARTINI/run_sim_water.sh` so Slurm submission is a single entry point.
- Use the existing draft scripts (`scan_water_diffusion.py`, `plot_diffusion_results.py`) and refine them rather than introducing new tooling.
- Compute diffusion only from the last 50% of the production trajectory to match analysis requirements.

## Execution Phases
- [x] Inspect current scripts and workflow inputs/outputs for MARTINI water scan.
- [x] Update `scan_water_diffusion.py` to run the parameter sweep and capture results.
- [x] Update `plot_diffusion_results.py` to read results and generate the two plots.
- [x] Align `run_sim_water.sh` to orchestrate full workflow end-to-end.
- [ ] Sanity check the workflow for expected outputs and paths.

## Known Errors / Blockers
- None yet.
