import sys, os, shutil
import subprocess as sp
import numpy as np
import tables as tb
from math import sqrt
from pathlib import Path

upside_path = Path(os.environ.get("UPSIDE_HOME", Path(__file__).resolve().parents[2])).resolve()
upside_utils_dir = upside_path / "py"
sys.path.insert(0, str(upside_utils_dir))
import run_upside as ru
from martini_prepare_system import stage_seed

#----------------------------------------------------------------------
## Inputs (edit these)
#----------------------------------------------------------------------

pdb_id             = "glpG-RKRK-79HIS"
protein_aa_pdb     = os.path.expanduser("~/Downloads/glpG-RKRK-79HIS.pdb")        # all-atom protein
lipid_name         = "POPE,POPG"                                                  # bilayer moleculetypes in dryMARTINI_itp/
charmm_gui_dir     = os.path.expanduser("~/Downloads/charmm-gui-8662225279")      # CHARMM-GUI Martini POPE/POPG membrane
opm_reference      = os.path.expanduser("~/Downloads/2nr9.pdb")                   # OPM membrane-oriented reference
membrane_thickness = 48.8   # equilibrated dry-MARTINI membrane thickness (A) for the ion count

sim_id             = "REMD"
work_dir           = "./"

#----------------------------------------------------------------------
## Stage 1 settings: build the equilibrated production seed
#----------------------------------------------------------------------

salt_molar         = 0.15
protein_pbc_margin = 15.0    # PBC-safe belt around the protein (>= dry-MARTINI nonbonded cutoff)
prep_temperature   = 0.8647  # equilibration temperature (T_up)
eq_time_step       = 0.009
eq_frame_steps     = 1000
eq_stage_nsteps    = 500     # steps per graded head-restraint-release stage (6.2 -> 6.5)
eq_66_nsteps       = 60000   # 6.6: the box converges in 30,000 steps and the lipid face in 60,000 (findings 8b)
min_70_max_iter    = 500
burnin_70_nsteps   = 40000   # NVT burn-in, protein rigid; its final state is the seed (stage 7.0 runs no production)
prep_seed          = 2026

#----------------------------------------------------------------------
## Stage 2 settings: temperature-REMD HDX production
#----------------------------------------------------------------------

duration           = 200000  # production length (upside time units); extend with continue_sim = True
frame_interval     = 100     # a frame every 100 time units (~2000 frames) -- the Upside standard (see
                             # 04.HDX). Momentum is NOT recorded and restarts are position-only, so the
                             # trajectories stay small (~120 MB/replica) and never blow up the disk.
prod_time_step     = 0.009
n_rep              = 48      # replicas
T_low              = 0.70    # ladder bottom (T_up); target 0.8647 sits near the top
T_high             = 0.90    # ladder top (T_up)
replica_interval   = 10      # time between exchange attempts (upside time units)
continue_sim       = False   # True -> extend from the last frame of the existing replica trajectories
randomseed         = 1

account            = "pi-trsosnic"
slurm_args_file    = "/project/trsosnic/yinhan/slurm/midway2.args"   # midway2 partition and node exclusions (CLAUDE.md)
cpus_per_node      = 28              # broadwl; the replicas share the node's cores as OpenMP threads
run_time           = "35:45:00"      # requested wall time (hh:mm:ss)
job_name           = "{}_{}".format(pdb_id, sim_id)

#----------------------------------------------------------------------
# Set the paths and filenames
#----------------------------------------------------------------------

script_dir = Path(__file__).resolve().parent
work_dir_path = Path(work_dir).expanduser()
if not work_dir_path.is_absolute():
    work_dir_path = script_dir / work_dir_path
prep_dir = work_dir_path / "outputs" / "martini_{}_hybrid".format(pdb_id)
run_dir  = prep_dir / sim_id
seed_up  = prep_dir / "checkpoints" / "{}.stage_7.0.up".format(pdb_id)
run_dir.mkdir(parents=True, exist_ok=True)

h5_files = [run_dir / "{}.run.{}.up".format(pdb_id, j) for j in range(n_rep)]
h5_files_str = " ".join(str(h5) for h5 in h5_files)
log_file = run_dir / "{}.run.log".format(pdb_id)

python = upside_path / ".venv" / "bin" / "python3"
if not python.exists():
    python = Path(sys.executable)

#----------------------------------------------------------------------
## Stage 1: prepare the equilibrated production seed
#----------------------------------------------------------------------
# CHARMM-GUI membrane + martinize CG envelope + OPM orientation + xy-only Monte-Carlo-barostat NPT +
# solvent-volume ions -> one equilibrated production .up. Runs once per variant (see readme.md).

if not continue_sim and not seed_up.exists():
    param_file = upside_path / "parameters" / "ff_2.1" / "martini.h5"
    if not param_file.exists():
        print("Generating the dry-MARTINI force-field file...")
        sp.check_call([str(python), str(upside_utils_dir / "martini_build_tables.py"),
                       "--upside-home", str(upside_path)])

    print("Preparing the hybrid seed...")
    prep_cmd = [
        str(python), str(upside_utils_dir / "martini_prepare_system.py"), "run-hybrid-workflow",
        "--pdb-id", pdb_id,
        "--upside-home", str(upside_path),
        "--run-dir", str(prep_dir),
        "--protein-aa-pdb", str(Path(protein_aa_pdb).expanduser().resolve()),
        "--lipid-name", lipid_name,
        "--charmm-gui-dir", str(Path(charmm_gui_dir).expanduser().resolve()),
        "--protein-orientation-mode", "opm",
        "--opm-reference", str(Path(opm_reference).expanduser().resolve()),
        "--membrane-thickness-angstrom", str(membrane_thickness),
        "--salt-molar", str(salt_molar),
        "--protein-pbc-margin", str(protein_pbc_margin),
        "--temperature", str(prep_temperature),
        "--min-70-max-iter", str(min_70_max_iter),
        "--eq-62-nsteps", str(eq_stage_nsteps),
        "--eq-63-nsteps", str(eq_stage_nsteps),
        "--eq-64-nsteps", str(eq_stage_nsteps),
        "--eq-65-nsteps", str(eq_stage_nsteps),
        "--eq-66-nsteps", str(eq_66_nsteps),
        "--prod-70-burnin-nsteps", str(burnin_70_nsteps),
        "--prod-70-nsteps", "0",
        "--prod-70-npt-enable", "0",
        "--eq-time-step", str(eq_time_step),
        "--prod-time-step", str(prod_time_step),
        "--eq-frame-steps", str(eq_frame_steps),
        "--prep-seed", str(prep_seed),
        "--seed", str(randomseed),
    ]
    print(" ".join(prep_cmd))
    sp.check_call(prep_cmd, cwd=str(script_dir))

if not seed_up.exists():
    raise FileNotFoundError("Production seed not found: {}\n"
                            "Run Stage 1 (continue_sim = False) first.".format(seed_up))

#----------------------------------------------------------------------
## Stage 2: replicate the seed (or continue) and launch temperature-REMD
#----------------------------------------------------------------------

# Each segment draws its own noise: a position-only restart (momentum is not recorded, to keep the files
# small) under the same seed would replay the previous segment's random numbers from round 0.
segment = 0
if continue_sim:
    print("Archiving previous output and reseeding from the last frame...")
    for fn in h5_files:
        if not fn.exists():
            raise FileNotFoundError("continue_sim = True but missing {}".format(fn))
        with tb.open_file(str(fn), "a") as t:
            i = 0
            while "output_previous_%i" % i in t.root:
                i += 1
            last = t.root.output if "output" in t.root else t.get_node("/output_previous_%i" % (i - 1))
            t.root.input.pos[:, :, 0] = last.pos[-1, 0]
            if "output" in t.root:
                t.root.output._f_rename("output_previous_%i" % i)
                i += 1
        segment = i
else:
    for fn in h5_files:
        shutil.copyfile(str(seed_up), str(fn))
segment_seed = stage_seed(randomseed, "remd.segment%i" % segment)

# Temperature ladder, uniform in sqrt(T) for near-constant exchange acceptance.
tempers = np.linspace(sqrt(T_low), sqrt(T_high), n_rep) ** 2
tempers_str = ",".join(str(t) for t in tempers)
swap_sets = ru.swap_table2d(1, n_rep)   # alternating even/odd neighbour pairs

upside_opts = (
    "--duration {} --frame-interval {} --temperature {} "
    "--time-step {} --integrator v --thermostat-timescale 5.0 --thermostat-interval -1 "
    "--seed {} --disable-recentering "
    "--replica-interval {} --swap-set {} --swap-set {} "
).format(duration, frame_interval, tempers_str, prod_time_step, segment_seed,
         replica_interval, swap_sets[0], swap_sets[1])

n_threads = min(n_rep, cpus_per_node)
sbatch_opts = (
    "$(cat {}) --account={} --job-name={} --output={} --time={} "
    "--nodes=1 --ntasks-per-node=1 --cpus-per-task={} "
).format(slurm_args_file, account, job_name, log_file, run_time, n_threads)

print("Launching REMD ({} replicas, T {}-{}, segment {})...".format(n_rep, T_low, T_high, segment))
cmd = 'sbatch {} --wrap="OMP_NUM_THREADS={} {}/obj/upside {} {}"'.format(
    sbatch_opts, n_threads, upside_path, upside_opts, h5_files_str)
print(cmd)
sp.check_call(cmd, shell=True)
