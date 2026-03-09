#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH -q batch
#SBATCH -t 08:00:00
#SBATCH -A marine-cpu
#SBATCH -J eval_stat
#SBATCH -o logs/eval_stat.out

set -euo pipefail

module use /work2/noaa/marine/jmeixner/general/modulefiles
module load ww3tools

WORKDIR="/work2/noaa/marine/ming.chen/GFS_Retro_Data"
WW3TOOLSDIR="${WORKDIR}/WW3-tools/ww3tools"

export PYTHONPATH="${WW3TOOLSDIR}:${PYTHONPATH:-}"

PY_SCRIPT="eval_stat.py"

python -u "${PY_SCRIPT}"
