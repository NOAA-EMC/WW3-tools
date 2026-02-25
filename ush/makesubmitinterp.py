import datetime as dt
from dateutil.relativedelta import relativedelta
import os
import re
import glob

## ===================== Setting (modified as needed) =========================
MACHINE = "ursa" # machine name ursa/orion/hercules
WORKDIR = "/scratch3/NCEPDEV/marine/Ming.Chen/ursa/ww3tools"

MODEL_BASE = "/scratch3/NCEPDEV/climate/Jessica.Meixner/Data/gfsv16"
SAT_BASE = "/scratch3/NCEPDEV/climate/Jessica.Meixner/WaveEvaluation/processsatdata/combineoutmonthly"   # if empty, the default directory will be used as WORKDIR/processsatdata/combineoutmonthly

# satellite and model settings
satellites=['JASON3', 'CRYOSAT2']
model="GFSv16" # now only support model of GFSv16 and retrov17_01
tz_list = ["00","06","12","18"]
grid = "global.0p25"

# Slurm settings
SBATCH_ACCOUNT   = "marine-cpu"
SBATCH_QUEUE     = "batch"
SBATCH_TIME      = "08:00:00"
SBATCH_NODES         = 1
SBATCH_NTASKS        = 1
SBATCH_CPUS_PER_TASK = 4
SBATCH_MEM           = "16G"
SET_THREAD_ENVS = True

## --------------- Machine-specific configuration ------------------

if MACHINE == "ursa":
    MODULE_USE_PATH = "/scratch4/NCEPDEV/marine/Saeideh.Banihashemi/installs/python-modules/"
    MODULE_LOAD     = "Ursa_ENV"
elif MACHINE in ("orion", "hercules"):
    MODULE_USE_PATH = "/work2/noaa/marine/jmeixner/general/modulefiles"
    MODULE_LOAD     = "ww3tools"
else:
    print(f"ERROR: Unsupported MACHINE='{MACHINE}'. Use Ursa, Orion, or Hercules.", file=sys.stderr)
    sys.exit(1)

## --------------- Directory settings ------------------------------
rootdir = os.path.join(WORKDIR, "processsatdata", "jobinterp")

if not SAT_BASE or not SAT_BASE.strip():
    SAT_BASE = os.path.join(WORKDIR, "processsatdata", "combineoutmonthly")

OUTDIR_BASE = os.path.join(WORKDIR, "processsatdata", "outinterp", model)

PROC_SCRIPT = os.path.join(WORKDIR, "WW3-tools", "ww3tools", "ProcSat_interpolation.py")

## --------------- Checking inputs and settings --------------------

jobdir = os.path.join(rootdir, model)

if not os.path.isdir(rootdir):
    os.makedirs(rootdir, exist_ok=True)

if not os.path.isdir(jobdir):
    os.makedirs(jobdir, exist_ok=True)

if not os.path.isdir(OUTDIR_BASE):
    os.makedirs(OUTDIR_BASE, exist_ok=True)

re_gfs = re.compile(r"^gfs\.(\d{8})$")
re_sat = re.compile(r"^Altimeter_(.+)_(\d{6})\.nc$")

def discover_model_dates(model_base: str):
    """
    Discover available model dates from folders:
      gfs.YYYYMMDD
    Returns a sorted list of YYYYMMDD strings.
    """
    dates = []
    for name in sorted(os.listdir(model_base)):
        if re.match(r"^gfs\.\d{8}$", name):
            yyyymmdd = name.split(".")[1]
            dates.append(yyyymmdd)
    return sorted(dates)

def sat_month_available_all(sat_base: str, satellites, yyyymm: str) -> bool:
    """
    Return True only if ALL satellites have monthly files for yyyymm:
      Altimeter_{sat}_{yyyymm}.nc
    """
    for sat in satellites:
        fname = os.path.join(sat_base, f"Altimeter_{sat}_{yyyymm}.nc")
        if not os.path.isfile(fname):
            return False
    return True


cdates = discover_model_dates(MODEL_BASE)

total = len(cdates)
covered = 0
missing = 0

for cdate in cdates:
    yyyymm = cdate[:6]
    if sat_month_available_all(SAT_BASE, satellites, yyyymm):
        covered += 1
    else:
        missing += 1

print("Satellite coverage verification:")
print(f"  Total model dates        : {total}")
print(f"  Dates fully covered      : {covered}")
print(f"  Dates missing satellites : {missing}")

# ------------------- Write jobcards -----------------------------------
written = 0
skipped_no_sat_all = 0
skipped_no_model_gribs = 0
missing_cycles = []

for cdate in cdates:
    yyyymm = cdate[:6]

    sat_ok = True
    for sat in satellites:
        fname = os.path.join(SAT_BASE, f"Altimeter_{sat}_{yyyymm}.nc")
        if not os.path.isfile(fname):
            sat_ok = False
            break
    if not sat_ok:
        skipped_no_sat_all += 1
        continue

    for tz in tz_list:
        if model == "GFSv16":
            model_gridded_dir = os.path.join(MODEL_BASE, f"gfs.{cdate}", tz, "wave", "gridded")
            MODEL_DATA_PATTERN_TEMPLATE = "gfswave.t{tz}z.{grid}.f*.grib2"
        elif model == "retrov17_01":
            model_gridded_dir = os.path.join(MODEL_BASE, f"gfs.{cdate}", tz, "products", "wave", "gridded","global.0p25")
            MODEL_DATA_PATTERN_TEMPLATE = "gfs.t{tz}z.{grid}.f*.grib2"
        else:
            print(f"ERROR: Unsupported Model.", file=sys.stderr)
            sys.exit(1)

        if not os.path.isdir(model_gridded_dir):
            skipped_no_model_gribs += 1
            missing_cycles.append(f"{cdate}{tz}")
            continue

        pattern = MODEL_DATA_PATTERN_TEMPLATE.format(tz=tz, grid=grid)
        gribs = glob.glob(os.path.join(model_gridded_dir, pattern))
        if len(gribs) == 0:
            skipped_no_model_gribs += 1
            missing_cycles.append(f"{cdate}{tz}")
            continue

        cdate_full = f"{cdate}{tz}"
        outfile = os.path.join(jobdir, f"job_{model}_{grid}_{cdate_full}.sh")
        outlog  = os.path.join(jobdir, f"run_{model}_{grid}_{cdate_full}.o%j")

        with open(outfile, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(f"#SBATCH --nodes={SBATCH_NODES}\n")
            f.write(f"#SBATCH --ntasks={SBATCH_NTASKS}\n")
            f.write(f"#SBATCH --cpus-per-task={SBATCH_CPUS_PER_TASK}\n")
            f.write(f"#SBATCH --mem={SBATCH_MEM}\n")
            f.write(f"#SBATCH -q {SBATCH_QUEUE}\n")
            f.write(f"#SBATCH -t {SBATCH_TIME}\n")
            f.write(f"#SBATCH -A {SBATCH_ACCOUNT}\n")
            f.write(f"#SBATCH -J procsat_{model}_{grid}_{cdate_full}\n")
            f.write(f"#SBATCH -o {outlog}\n")

            f.write(f"module use {MODULE_USE_PATH}\n")
            f.write(f"module load {MODULE_LOAD}\n\n")

            f.write("set -euo pipefail\n")
            f.write("set -x\n\n")

            if SET_THREAD_ENVS:
                f.write("export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n")
                f.write("export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK\n\n")

            f.write(f"MODEL={model}\n")
            f.write(f"GRID={grid}\n")
            f.write(f"DATE={cdate}\n")
            f.write(f"TZ={tz}\n")
            f.write(f"YYYYMM={yyyymm}\n")
            f.write(f"OUTDIR={OUTDIR_BASE}\n")
            f.write("mkdir -p ${OUTDIR}\n\n")

            if model == "GFSv16":
                f.write(f"MODEL_DATA_DIR={MODEL_BASE}/gfs.${{DATE}}/${{TZ}}/wave/gridded\n")
            elif model == "retrov17_01":
                f.write(f"MODEL_DATA_DIR={MODEL_BASE}/gfs.${{DATE}}/${{TZ}}/products/wave/gridded/global.0p25\n")
            else:
                print(f"ERROR: Unsupported Model.", file=sys.stderr)
                sys.exit(1)

            f.write(f"MODEL_DATA_PATTERN='{pattern}'\n\n")

            for sat in satellites:
                f.write(f"SAT={sat}\n")
                f.write(f"SATELLITE_FILE={SAT_BASE}/Altimeter_{sat}_${{YYYYMM}}.nc\n")
                f.write(f"OUTPUT_FILE=${{MODEL}}_${{GRID}}_{cdate_full}_{sat}.nc\n")
                f.write(
                    f"python {PROC_SCRIPT} "
                    f"-t grib2 -d $MODEL_DATA_DIR -p $MODEL_DATA_PATTERN "
                    f"-s $SATELLITE_FILE -o $OUTDIR -f $OUTPUT_FILE -m $MODEL\n\n"
                )

        os.chmod(outfile, 0o750)
        written += 1

print("\nJobcard generation summary:")
print(f"  Jobcards written                              : {written}")
print(f"  Dates skipped (month missing ≥1 satellite)     : {skipped_no_sat_all}")
print(f"  Cycles skipped (missing model dir or GRIBs)    : {skipped_no_model_gribs}")

if missing_cycles:
    for c in sorted(missing_cycles):
        print(f"    {c}")

ush = os.path.join(WORKDIR, "WW3-tools", "ush", "run_all_jobs.sh")
os.makedirs(jobdir, exist_ok=True)
cmd = f"cp {ush} {jobdir}"
os.system(cmd)


print(f"  Jobcards directory                             : {jobdir}")
