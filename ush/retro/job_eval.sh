#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH -q batch
#SBATCH -t 08:00:00
#SBATCH -A marine-cpu
#SBATCH -J eval
#SBATCH -o logs/eval.out

set -euo pipefail

module use /work2/noaa/marine/jmeixner/general/modulefiles
module load ww3tools

WORKDIR="/work2/noaa/marine/ming.chen/GFS_Retro_Data"
WW3TOOLSDIR="${WORKDIR}/WW3-Tools/ww3tools"
CURRENTDIR=$(pwd)

FILES=("mvalstats.py" "pvalstats.py")

FILES=("mvalstats.py" "pvalstats.py")

for f in "${FILES[@]}"; do
    if [[ -f "${CURRENTDIR}/${f}" ]]; then
        echo "${f} exists in current directory."
    else
        echo "${f} not found. Copying from ${WW3TOOLSDIR}..."
        cp "${WW3TOOLSDIR}/${f}" "${CURRENTDIR}/"
    fi
done

PY_SCRIPT="eval.py"

python -u "${PY_SCRIPT}"
