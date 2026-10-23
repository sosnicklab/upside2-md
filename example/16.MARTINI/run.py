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
from martini_prepare_system import stage_seed, promote_md_output_state_to_input, BELT_GATE_PRODUCTION_STEPS

#----------------------------------------------------------------------
## Inputs (edit these)
#----------------------------------------------------------------------

pdb_id             = "glpG-RKRK-79HIS"
protein_aa_pdb     = os.path.expanduser("~/Downloads/glpG-RKRK-79HIS.pdb")        # all-atom protein
lipid_name         = "POPE,POPG"                                                  # bilayer moleculetypes in dryMARTINI_itp/
charmm_gui_dir     = os.path.expanduser("~/Downloads/charmm-gui-8662225279")      # CHARMM-GUI Martini POPE/POPG membrane
opm_reference      = os.path.expanduser("~/Downloads/2nr9.pdb")                   # OPM membrane-oriented reference
membrane_thickness = 41.8   # POPE/POPG P-P thickness (A) of the GROMACS reference bilayer (findings 8), for the ion count

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
burnin_70_nsteps   = 40000   # NVT burn-in, protein rigid; then stage 7.0's production carries the solvation gate
prep_seed          = 2026

#----------------------------------------------------------------------
## Stage 2 settings: temperature-REMD HDX production
#----------------------------------------------------------------------

n_blocks           = 4       # REMD blocks of one wall time each; raise it and re-run to extend
frame_interval     = 10      # time units per frame: a block of 48 glpG replicas on 28 cores covers ~1,000 time units
                             # (~28 per wall hour), so ~100 frames per block. Momentum is not recorded and blocks
                             # restart from positions, so a replica grows by one frame of positions per frame.
prod_time_step     = 0.009
n_rep              = 48      # replicas
T_low              = 0.70    # ladder bottom (T_up); target 0.8647 sits near the top
T_high             = 0.90    # ladder top (T_up)
replica_interval   = 10      # time between exchange attempts (upside time units)
randomseed         = 1

account            = "pi-trsosnic"
slurm_args_file    = "/project/trsosnic/yinhan/slurm/midway2.args"   # midway2 partition and node exclusions (CLAUDE.md)
job_env_file       = "/beagle3/trsosnic/yinhan/upside3/env_upside3.sh"   # each job sources the tree's environment
cpus_per_node      = 28              # broadwl; the replicas share the node's cores as OpenMP threads
run_time           = "35:45:00"      # requested wall time of every job (hh:mm:ss)
time_margin        = 1800            # seconds a block keeps for its final frame, archiving and the next submission
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
block_count = run_dir / "block_count"
n_threads = min(n_rep, cpus_per_node)

# A job runs this script again with HYBRID_JOB_STAGE naming its stage; run by hand, it submits the next job.
job_stage = os.environ.get("HYBRID_JOB_STAGE", "")


def submit(stage, log, dependency=None):
    cmd = ('sbatch --parsable $(cat {}) --account={} --job-name={}_{} --output={} --time={} '
           '--nodes=1 --ntasks-per-node=1 --cpus-per-task={} {}--wrap=". {} && export OMP_NUM_THREADS={} '
           'HYBRID_JOB_STAGE={} && cd {} && python3 {}"').format(
        slurm_args_file, account, job_name, stage, log, run_time, n_threads,
        '--dependency=afterok:{} '.format(dependency) if dependency else '',
        job_env_file, n_threads, stage, script_dir, Path(__file__).resolve())
    job = sp.check_output(cmd, shell=True, text=True).strip()
    print("submitted {} job {}".format(stage, job))
    return job

#----------------------------------------------------------------------
## Stage 1: prepare the equilibrated production seed
#----------------------------------------------------------------------
# CHARMM-GUI membrane + martinize CG envelope + OPM orientation + xy-only Monte-Carlo-barostat NPT +
# solvent-volume ions -> one equilibrated production .up. Runs once per variant (see readme.md).

if job_stage == "prep":
    param_file = upside_path / "parameters" / "ff_2.1" / "martini.h5"
    if not param_file.exists():
        print("Generating the dry-MARTINI force-field file...")
        sp.check_call([sys.executable, str(upside_utils_dir / "martini_build_tables.py"),
                       "--upside-home", str(upside_path)])

    print("Preparing the hybrid seed...")
    prep_cmd = [
        sys.executable, str(upside_utils_dir / "martini_prepare_system.py"), "run-hybrid-workflow",
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
        "--prod-70-nsteps", str(BELT_GATE_PRODUCTION_STEPS),
        "--prod-70-npt-enable", "0",
        "--eq-time-step", str(eq_time_step),
        "--prod-time-step", str(prod_time_step),
        "--eq-frame-steps", str(eq_frame_steps),
        "--prep-seed", str(prep_seed),
        "--seed", str(randomseed),
    ]
    print(" ".join(prep_cmd))
    sp.check_call(prep_cmd, cwd=str(script_dir))
    sys.exit(0)

#----------------------------------------------------------------------
## Stage 2: one block of temperature-REMD
#----------------------------------------------------------------------
# The first block replicates the seed, whose stage-7.0 production passed the solvation gate, and starts every
# replica from its last frame; each later block restarts each slot from its own last frame with the walker label it
# held, its output archived as output_previous_N. Each block draws its own noise: momentum is not recorded, so a
# position-only restart under the same seed would replay the previous block's random numbers from round 0. The
# engine stops itself at --time-limit, before the wall, and the block queues the next one up to n_blocks.

if job_stage == "remd":
    block = int(block_count.read_text()) + 1 if block_count.exists() else 1   # block_count counts completed blocks
    segment = 0
    if block > 1:
        for fn in h5_files:
            with tb.open_file(str(fn), "a") as t:
                i = 0
                while "output_previous_%i" % i in t.root:
                    i += 1
                last = t.root.output if "output" in t.root else t.get_node("/output_previous_%i" % (i - 1))
                t.root.input.pos[:, :, 0] = last.pos[-1, 0]
                t.root.input.mom._v_attrs.replica_index = int(last.replica_index[-1, 0])  # the walker the slot held
                if "output" in t.root:
                    t.root.output._f_rename("output_previous_%i" % i)
                    i += 1
            segment = i
    else:
        start_up = run_dir / "{}.remd_start.up".format(pdb_id)
        shutil.copyfile(str(seed_up), str(start_up))
        promote_md_output_state_to_input(start_up)
        for fn in h5_files:
            shutil.copyfile(str(start_up), str(fn))

    # Temperature ladder, uniform in sqrt(T) for near-constant exchange acceptance.
    tempers = np.linspace(sqrt(T_low), sqrt(T_high), n_rep) ** 2
    swap_sets = ru.swap_table2d(1, n_rep)   # alternating even/odd neighbour pairs (one set for two replicas)
    h, m, s = (int(x) for x in run_time.split(":"))
    upside_cmd = [str(upside_path / "obj" / "upside")] + [str(fn) for fn in h5_files] + [
        "--duration", "1e6", "--time-limit", str(3600 * h + 60 * m + s - time_margin),
        "--frame-interval", str(frame_interval), "--temperature", ",".join(str(t) for t in tempers),
        "--time-step", str(prod_time_step), "--integrator", "v", "--thermostat-timescale", "5.0",
        "--thermostat-interval", "-1", "--seed", str(stage_seed(randomseed, "remd.segment%i" % segment)),
        "--disable-recentering", "--replica-interval", str(replica_interval)] + [
        x for swap_set in swap_sets for x in ("--swap-set", swap_set)]
    print("REMD block {} of {} ({} replicas, T {}-{}, segment {})".format(block, n_blocks, n_rep, T_low, T_high, segment))
    sp.check_call(upside_cmd)
    block_count.write_text(str(block))
    if block < n_blocks:
        submit("remd", run_dir / "{}.remd_block{}.log".format(pdb_id, block + 1))
    sys.exit(0)

#----------------------------------------------------------------------
## By hand: submit what comes next
#----------------------------------------------------------------------

if not seed_up.exists():
    prep_job = submit("prep", prep_dir / "{}.prep.log".format(pdb_id))
    submit("remd", run_dir / "{}.remd_block1.log".format(pdb_id), dependency=prep_job)
else:
    done = int(block_count.read_text()) if block_count.exists() else 0
    if done >= n_blocks:
        print("{} of {} REMD blocks done; raise n_blocks to extend".format(done, n_blocks))
    else:
        submit("remd", run_dir / "{}.remd_block{}.log".format(pdb_id, done + 1))
