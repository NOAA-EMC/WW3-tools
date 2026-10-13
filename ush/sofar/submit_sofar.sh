#!/bin/sh
#SBATCH --nodes=1
#SBATCH -q debug
#SBATCH -t 00:30:00
#SBATCH -A marine-cpu
#SBATCH -J submit_sofar
#SBATCH -o submit_sofar.out
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32        # Adjust as the parallelization is improved
##SBATCH --mem=0                  # Request all available memory on the node
##SBATCH --exclusive
#SBATCH --open-mode=truncate

# ================= USER'S INPUT =======================
START_DATE="2026-01-01"
END_DATE="2026-01-02"
WORK_DIR="/scratch3/NCEPDEV/marine/Edwin.Rainville/ww3-tools-dev/data"
ENV_PATH="/scratch3/NCEPDEV/marine/Edwin.Rainville/envs/sofar"
# ======================================================

# Load RDHPCS-conda and activate environment
module load rdhpcs-conda
eval "$(conda shell.bash hook)"
conda activate $ENV_PATH

# Print out the user settings to the log 
echo "==== Submit SOFAR Validation Job ===="
SUBMIT_TIME=$(scontrol show job $SLURM_JOB_ID | grep -o 'SubmitTime=[^ ]*' | cut -d= -f2)
echo "Job Submitted: $SUBMIT_TIME"
echo "Validation Start Date: $START_DATE"
echo "Validation End Date: $END_DATE"
echo "Work Directory: $WORK_DIR"
echo "Environment path: $ENV_PATH"

# Execute the python matchup script
python sofar_matchup_pipeline.py \
    --start-date "${START_DATE}" \
    --end-date "${END_DATE}" \
    --work-dir "${WORK_DIR}"
