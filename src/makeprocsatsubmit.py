import datetime as dt
from dateutil.relativedelta import relativedelta
import os
import sys

# ================================================
#           USER-EDITABLE CONFIGURATION
# ================================================
MACHINE = "ursa" # or orion/herculus
WORKDIR = "/scratch3/NCEPDEV/marine/Ming.Chen/ursa/ww3tools" # working directory including WW3-tools and processeddata

STARTDATE = "2024-11-15" # start date with formats YYYY-MM-DD or YYYYMMDD
ENDDATE   = "2025-01-15" # end date with formats YYYY-MM-DD or YYYYMMDD

SATELLITES = "JASON3,CRYOSAT2,SARAL,SENTINEL3A" # satellites using comma or space separated

# SLURM settings
SBATCH_QUEUE     = "batch"
SBATCH_ACCOUNT   = "marine-cpu"
SBATCH_WALLTIME  = "08:00:00"

SBATCH_EXCLUSIVE = True # True: exclusive mode (whole node)

  # if SBATCH_EXCLUSIVE = True, the settings below are ignored
SBATCH_NODES         = "1"
SBATCH_NTASKS         = "1"
SBATCH_CPUS_PER_TASK  = "4"
SBATCH_MEM            = "16G"

SET_THREAD_ENVS       = True # Set OMP_NUM_THREADS, MKL_NUM_THREADS, etc.

# ===============================================

ROOTDIR = os.path.join(WORKDIR, "processsatdata", "jobsubs")    # output jobcards directory
THISDIR = os.path.join(WORKDIR, "WW3-tools", "src")             # source script directory
PATHTOWW3TOOLS = os.path.join(WORKDIR, "WW3-tools", "ww3tools") # ww3tools directory (ProcSat_Altimeter.py)
OUT_BASE = os.path.join(WORKDIR, "processsatdata", "out")       # output directory for processed data (origional defined in .yaml)

MACHINE = MACHINE.strip().lower()

# Module commands (different machines have different module builds and directories)
if MACHINE == "ursa":
    MODULE_USE_PATH = "/scratch3/NCEPDEV/climate/Jessica.Meixner/general/modulefiles"
    MODULE_LOAD     = "ww3tools"
    YAML_CONFIG_SUBDIR = "configs_ursa"

elif MACHINE in ("orion", "hercules"):
    MODULE_USE_PATH = "/work2/noaa/marine/jmeixner/general/modulefiles"
    MODULE_LOAD     = "ww3tools"
    YAML_CONFIG_SUBDIR = "configs_orion"

else:
    print(f"ERROR: Unsupported MACHINE='{MACHINE}'. Use Ursa, Orion, or Hercules.", file=sys.stderr)
    sys.exit(1)

# Check required directories exist
if not os.path.isdir(THISDIR):
    print(f"ERROR: THISDIR does not exist: {THISDIR}", file=sys.stderr)
    sys.exit(1)

if not os.path.isdir(PATHTOWW3TOOLS):
    print(f"ERROR: PATHTOWW3TOOLS does not exist: {PATHTOWW3TOOLS}", file=sys.stderr)
    sys.exit(1)

# Create output directory if missing
if not os.path.isdir(ROOTDIR):
    try:
        os.makedirs(ROOTDIR)
        print(f"Created output directory: {ROOTDIR}")
    except Exception as e:
        print(f"ERROR: Cannot create ROOTDIR: {e}", file=sys.stderr)
        sys.exit(1)

# Parse satellite list
sats = [s.strip() for s in SATELLITES.replace(",", " ").split() if s.strip()]
if not sats:
    print("ERROR: No valid satellites provided", file=sys.stderr)
    sys.exit(1)

# Parse dates
def parse_date(s):
    s = s.strip()
    for fmt in ("%Y-%m-%d", "%Y%m%d"):
        try:
            return dt.datetime.strptime(s, fmt)
        except ValueError:
            continue
    print(f"ERROR: Invalid date format: '{s}'. Use YYYY-MM-DD or YYYYMMDD.", file=sys.stderr)
    sys.exit(1)

# Generate date pairs (15-day steps + monthly overlap pattern)
dates1 = []
dates2 = []

current = parse_date(STARTDATE)
enddate = parse_date(ENDDATE)

while current <= enddate:
    d1 = current.strftime("%Y%m%d")
    d2 = (current + dt.timedelta(days=15)).strftime("%Y%m%d")
    dates1.append(d1)
    dates2.append(d2)

    dates1.append(d2)
    current += relativedelta(months=+1)
    dates2.append(current.strftime("%Y%m%d"))

job_count = 0

for i in range(len(dates1)):
    for sat in sats:
        jobname = f"job_{sat}_{dates1[i]}.sh"
        filepath = os.path.join(ROOTDIR, jobname)

        with open(filepath, "w", encoding="utf-8") as f:
            f.write("#!/bin/bash\n\n")

            # Common SLURM directives
            f.write(f"#SBATCH --nodes={SBATCH_NODES}\n")
            f.write(f"#SBATCH -q {SBATCH_QUEUE}\n")
            f.write(f"#SBATCH -t {SBATCH_WALLTIME}\n")
            f.write(f"#SBATCH -A {SBATCH_ACCOUNT}\n")
            f.write(f"#SBATCH -J procsat_{sat}_{dates1[i]}\n")
            f.write(f"#SBATCH -o run_{sat}_{dates1[i]}.o%j\n")

            # Exclusive vs explicit resources
            if SBATCH_EXCLUSIVE:
                f.write("#SBATCH --exclusive\n")
            else:
                f.write(f"#SBATCH --ntasks={SBATCH_NTASKS}\n")
                f.write(f"#SBATCH --cpus-per-task={SBATCH_CPUS_PER_TASK}\n")
                f.write(f"#SBATCH --mem={SBATCH_MEM}\n")

            f.write("\n")

            # Module environment
            if MODULE_USE_PATH.strip():
                f.write(f"module use {MODULE_USE_PATH}\n")
            if MODULE_LOAD.strip():
                f.write(f"module load {MODULE_LOAD}\n")
            if MODULE_USE_PATH.strip() or MODULE_LOAD.strip():
                f.write("\n")

            # Thread environment control (only in shared mode)
            if SET_THREAD_ENVS and not SBATCH_EXCLUSIVE:
                f.write("# Control number of threads for performance & memory\n")
                f.write(f"export OMP_NUM_THREADS={SBATCH_CPUS_PER_TASK}\n")
                f.write(f"export MKL_NUM_THREADS={SBATCH_CPUS_PER_TASK}\n")
                f.write(f"export NUMEXPR_NUM_THREADS={SBATCH_CPUS_PER_TASK}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={SBATCH_CPUS_PER_TASK}\n")
                f.write("\n")

            # Job variables
            f.write(f'ThisDir="{THISDIR}"\n')
            f.write(f'PathToWW3TOOLS="{PATHTOWW3TOOLS}"\n')
            f.write(f'SAT="{sat}"\n')
            f.write(f'IDATE="{dates1[i]}00"\n')
            f.write(f'EDATE="{dates2[i]}00"\n\n')

            # The processing command
            f.write(f'YAMLFILE="${{ThisDir}}/{YAML_CONFIG_SUBDIR}/${{SAT}}.yaml"\n')
            f.write(f'OUT_BASE="{OUT_BASE}"\n\n')
            f.write('python "${PathToWW3TOOLS}/ProcSat_Altimeter.py" \\\n')
            f.write('    --satelite "${SAT}" \\\n')
            f.write('    --initdate "${IDATE}" \\\n')
            f.write('    --enddate "${EDATE}" \\\n')
            f.write('    --timestep 1.0 \\\n')
            f.write('    --yaml "${YAMLFILE}" \\\n')
            f.write('    --out_base "${OUT_BASE}"\n')

        # Make executable
        os.chmod(filepath, 0o755)
        job_count += 1
        print(f"Created: {jobname}")

ush = os.path.join(WORKDIR, "WW3-tools", "ush", "run_all_jobs.sh")
os.makedirs(ROOTDIR, exist_ok=True)
cmd = f"cp {ush} {ROOTDIR}"
os.system(cmd)

print(f"\nDone. Generated {job_count} job script(s) in:")
print(ROOTDIR)
