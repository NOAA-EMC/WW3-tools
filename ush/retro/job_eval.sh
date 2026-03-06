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

PY_SCRIPT="eval.py"

python -u "${PY_SCRIPT}"
