"""
makesubmitcombine.py

PURPOSE:
    Create Slurm job scripts to run `combineSatInterpOut.py`.
    This allows the combination tasks to run in parallel, improving
efficiency and reducing memory pressure compared with running all
model/satellite combinations in a single job.

USAGE:
 Modify the settings below as needed:
  - MACHINE:          HPC machine where the jobs will run.
                      Supported options: "ursa", "orion", "hercules"
  - WORKDIR:          Working directory containing `WW3-tools` and processed data.
  - MODELS:           List of model names to process.
                      These should be consistent with the model names used in the
                      interpolated output filenames.
  - SATELLITES:       List of satellite names to process.
  - STARTDATE:        Start date for the combination period, in `yyyymmddhh` format.
  - ENDDATE:          End date for the combination period, in `yyyymmddhh` format.
  - INTERVAL_HOURS:   Time interval, in hours, between forecast cycles to process.
                      for example: 00z 06z 12z 18z -> INTERVAL_HOURS=6
                                   00z 12z         -> INTERVAL_HOURS=12
                                   00z             -> INTERVAL_HOURS=24
  - MAX_FORECAST_DAY: Maximum forecast day to include in the combination.
  - FORCE_SEASON:     Label used in the output filename, such as a season, month,
                      or custom period name.
  - Slurm setting:    Define account, queue, wall-clock time, memory, CPUs, and
                      other job submission options.

OUTPUT:
 One Slurm job script is created for each model and satellite combination

 Output locaiton:
    {WORKDIR}/processsatdata/jobcombine/{model}/

 Output filename format:
     job_{model}_{satellite}_{force_season}.sh

AUTHOR and DATE:
 03/12/2026: Ming Chen, first version

"""

import datetime as dt
from dateutil.relativedelta import relativedelta
import os
import re
import glob
import sys
import shutil

## ===================== Setting (modified as needed) =========================
MACHINE = "orion" # machine name ursa/orion/hercules
WORKDIR = "/work2/noaa/marine/ming.chen/issue109"
MODELS  = ["GFSv16", "retrov17_01"]
#SATELLITES = ["JASON3", "CRYOSAT2", "SARAL", "SENTINEL3A", "SENTINEL3B", "SENTINEL6A"]
SATELLITES = ["JASON3"]
STARTDATE  = "2024120100"
ENDDATE    = "2025022800"

INTERVAL_HOURS   = 12
MAX_FORECAST_DAY = 16

FORCE_SEASON     = "DJF2025"

INPUTDIR_BASE = "/work2/noaa/marine/ming.chen/GFS_Retro_Data/data/outinterp"

# Slurm settings
SBATCH_ACCOUNT   = "marine-cpu"
SBATCH_QUEUE     = "batch"
SBATCH_TIME      = "08:00:00"
SBATCH_NODES         = 1
SBATCH_NTASKS        = 1
SBATCH_CPUS_PER_TASK = 4
SBATCH_MEM           = "180G"
SET_THREAD_ENVS = True

## --------------------- Machine-specific configuration -----------------------
if MACHINE == "ursa":
    MODULE_USE_PATH = "/scratch4/NCEPDEV/marine/Saeideh.Banihashemi/installs/python-modules/"
    MODULE_LOAD     = "Ursa_ENV"
elif MACHINE in ("orion", "hercules"):
    MODULE_USE_PATH = "/work2/noaa/marine/jmeixner/general/modulefiles"
    MODULE_LOAD     = "ww3tools"
else:
    print(f"ERROR: Unsupported MACHINE='{MACHINE}'. Use Ursa, Orion, or Hercules.", file=sys.stderr)
    sys.exit(1)

ROOTDIR = os.path.join(WORKDIR, "processsatdata", "jobcombine")
RUN_ALL_JOBS = os.path.join(WORKDIR, "WW3-tools", "ush", "run_all_jobs.sh")
SCRIPT_DIR = os.path.join(WORKDIR, "WW3-tools", "ush")
PY_SCRIPT = "combineSatInterpOut.py"

if not os.path.isfile(os.path.join(SCRIPT_DIR, PY_SCRIPT)):
    print(f"Error: Script not found: {SCRIPT_DIR}/{PY_SCRIPT}", file=sys.stderr)
    sys.exit(1)

def write_jobcard(model: str, satellite: str, outdir: str) -> str:
    """
    Write one combine jobcard for one model + one satellite.
    Returns the output filename.
    """
    # one jobcard per satellite
    outfile = os.path.join(outdir, f"job_{model}_{satellite}_{FORCE_SEASON}.sh")
    outlog = os.path.join(outdir, f"run_{model}_{satellite}_{FORCE_SEASON}.o%j")

    with open(outfile, "w") as f:
        f.write("#!/bin/bash\n")
        f.write(f"#SBATCH --nodes={SBATCH_NODES}\n")
        f.write(f"#SBATCH --ntasks={SBATCH_NTASKS}\n")
        f.write(f"#SBATCH --cpus-per-task={SBATCH_CPUS_PER_TASK}\n")
        f.write(f"#SBATCH --mem={SBATCH_MEM}\n")
        f.write(f"#SBATCH -q {SBATCH_QUEUE}\n")
        f.write(f"#SBATCH -t {SBATCH_TIME}\n")
        f.write(f"#SBATCH -A {SBATCH_ACCOUNT}\n")
        f.write(f"#SBATCH -J combine_{model}_{satellite}\n")
        f.write(f"#SBATCH -o {outlog}\n\n")

        f.write("set -euo pipefail\n")
        f.write("set -x\n\n")

        f.write("# =====================================================\n")
        f.write("# Load environment\n")
        f.write("# =====================================================\n\n")
        f.write(f"module use {MODULE_USE_PATH}\n")
        f.write(f"module load {MODULE_LOAD}\n\n")

        if SET_THREAD_ENVS:
            f.write("export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n")
            f.write("export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK\n\n")

        f.write("# =====================================================\n")
        f.write("# User Configuration\n")
        f.write("# =====================================================\n\n")
        f.write(f'WORKDIR="{WORKDIR}"\n')
        f.write(f'SCRIPT_DIR="{SCRIPT_DIR}"\n')
        f.write(f'PY_SCRIPT="{PY_SCRIPT}"\n\n')

        if INPUTDIR_BASE:
            f.write(f'INPUTDIR_BASE="{INPUTDIR_BASE}"\n\n')

        f.write(f'MODELS="{model}"\n')
        f.write(f'SATELLITES="{satellite}"\n\n')

        f.write(f'STARTDATE="{STARTDATE}"\n')
        f.write(f'ENDDATE="{ENDDATE}"\n\n')

        f.write(f"INTERVAL_HOURS={INTERVAL_HOURS}\n")
        f.write(f"MAX_FORECAST_DAY={MAX_FORECAST_DAY}\n\n")

        f.write(f'FORCE_SEASON="{FORCE_SEASON}"\n')

        f.write("# =====================================================\n")
        f.write("# Run\n")
        f.write("# =====================================================\n\n")

        f.write('echo "Running combineSatInterpOut.py"\n')
        f.write('echo "Models: ${MODELS}"\n')
        f.write('echo "Satellites: ${SATELLITES}"\n')
        f.write('echo "Start: ${STARTDATE}"\n')
        f.write('echo "End: ${ENDDATE}"\n')
        f.write('echo "WORKDIR: ${WORKDIR}"\n\n')

        # Keep downstream CLI compatible with your existing sample jobcard:
        f.write('python -u "${SCRIPT_DIR}/${PY_SCRIPT}" \\\n')
        f.write('    -models ${MODELS} \\\n')
        f.write('    -WORKDIR "${WORKDIR}" \\\n')
        f.write('    -satellites ${SATELLITES} \\\n')
        f.write('    -startdate "${STARTDATE}" \\\n')
        f.write('    -enddate "${ENDDATE}" \\\n')
        f.write('    -interval_hours ${INTERVAL_HOURS} \\\n')
        f.write('    -max_forecast_day ${MAX_FORECAST_DAY} \\\n')
        f.write('    -force_season ${FORCE_SEASON} \\\n')

        if INPUTDIR_BASE:
            f.write('    -INPUTDIR_BASE "${INPUTDIR_BASE}" \\\n')

    os.chmod(outfile, 0o750)
    return outfile

def main():

    written = 0
    if not os.path.isdir(ROOTDIR):
        os.makedirs(ROOTDIR, exist_ok=True)

    for model in MODELS:
        OUTDIR = os.path.join(ROOTDIR, model)
        os.makedirs(OUTDIR, exist_ok=True)

        for sat in SATELLITES:
            jobfile = write_jobcard(model, sat, OUTDIR)
            written += 1
            print(f"Wrote: {jobfile}")

        # copy run_all_jobs.sh into each model job directory
        if os.path.isfile(RUN_ALL_JOBS):
            dst = os.path.join(OUTDIR, os.path.basename(RUN_ALL_JOBS))
            shutil.copy2(RUN_ALL_JOBS, dst)
            print(f"Copied: {dst}")
        else:
            print(f"WARNING: run_all_jobs.sh not found: {RUN_ALL_JOBS}")

    print("\nJobcard generation summary:")
    print(f"  Models                : {len(MODELS)}")
    print(f"  Satellites per model  : {len(SATELLITES)}")
    print(f"  Total jobcards written: {written}")
    print(f"  Root output directory : {ROOTDIR}")


if __name__ == "__main__":
    main()
