import datetime as dt
from dateutil.relativedelta import relativedelta
import os
import sys

REQUIRED_VARS = [
    "ROOTDIR",
    "STARTDATE",
    "ENDDATE",
    "SATELITES",
    "THISDIR",
    "PATHTOWW3TOOLS",
]

def die(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)

def getenv_strict(name: str) -> str:
    val = os.environ.get(name, "").strip()
    if not val:
        die(f"Environment variable {name} is required but not set.")
    return val

def parse_date(s: str) -> dt.datetime:
    s = s.strip()
    # Accept YYYY-MM-DD or YYYYMMDD
    for fmt in ("%Y-%m-%d", "%Y%m%d"):
        try:
            return dt.datetime.strptime(s, fmt)
        except ValueError:
            pass
    die(f"Invalid date '{s}'. Use YYYY-MM-DD or YYYYMMDD.")

# ---- read all required variables strictly from env ----
for v in REQUIRED_VARS:
    if not os.environ.get(v):
        die(f"Environment variable {v} is required but not set.")

rootdir        = getenv_strict("ROOTDIR")
startdate_str  = getenv_strict("STARTDATE")
enddate_str    = getenv_strict("ENDDATE")
sat_env        = getenv_strict("SATELITES")
ThisDir        = getenv_strict("THISDIR")
PathToWW3TOOLS = getenv_strict("PATHTOWW3TOOLS")

startdate = parse_date(startdate_str)
enddate   = parse_date(enddate_str)

# satellites can be comma or space separated
satelites = [s for s in sat_env.replace(",", " ").split() if s]
if not satelites:
    die("SATELITES parsed to an empty list. Provide at least one satellite.")

nowdate = startdate
dates1 = []
dates2 = []

while nowdate <= enddate:
    dates1.append(nowdate.strftime('%Y%m%d'))
    dates2.append((nowdate + dt.timedelta(days=15)).strftime('%Y%m%d'))
    dates1.append((nowdate + dt.timedelta(days=15)).strftime('%Y%m%d'))
    nowdate = nowdate + relativedelta(months=+1)
    dates2.append(nowdate.strftime('%Y%m%d'))

# Generate job scripts (NO module lines inside, per your request)
for i in range(len(dates1)):
    for j in range(len(satelites)):
        jobname  = f"job_{satelites[j]}_{dates1[i]}.sh"
        outfile  = os.path.join(rootdir, jobname)
        with open(outfile, "w") as f:
            f.write("#!/bin/bash\n")
            f.write("#SBATCH --nodes=1\n")
            f.write("#SBATCH -q batch\n")
            f.write("#SBATCH -t 08:00:00\n")
            f.write("#SBATCH -A marine-cpu\n")
            f.write(f"#SBATCH -J procsat_{satelites[j]}_{dates1[i]}\n")
            f.write(f"#SBATCH -o run_{satelites[j]}_{dates1[i]}.o%j\n")
            f.write("#SBATCH --exclusive\n\n")
            f.write(f"module use /scratch4/NCEPDEV/marine/Saeideh.Banihashemi/installs/python-modules\n")
            f.write(f"module load Ursa_ENV\n\n")
            f.write(f'ThisDir="{ThisDir}"\n')
            f.write(f'PathToWW3TOOLS="{PathToWW3TOOLS}"\n\n')
            f.write(f'SAT="{satelites[j]}"\n')
            f.write(f'IDATE="{dates1[i]}00"\n')
            f.write(f'EDATE="{dates2[i]}00"\n\n')
            f.write('YAMLFILE="${ThisDir}/configs/${SAT}.yaml"\n')
            f.write('python "${PathToWW3TOOLS}/ProcSat_Altimeter.py" '
                    '--satelite "${SAT}" --initdate "${IDATE}" --enddate "${EDATE}" '
                    '--timestep 1.0 --yaml "${YAMLFILE}"\n')
        os.chmod(outfile, 0o755)
