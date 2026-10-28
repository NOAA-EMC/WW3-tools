#!/bin/sh
#SBATCH --nodes=1
#SBATCH -q debug
#SBATCH -t 00:30:00
#SBATCH -A marine-cpu
#SBATCH -J submit_sofar
#SBATCH -o submit_sofar.out
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --exclusive
#SBATCH --open-mode=truncate

# ================= USER'S INPUT =======================
START_DATE="2026-01-01"
END_DATE="2026-01-31"
WORK_DIR="/scratch3/NCEPDEV/marine/Edwin.Rainville/sofar_work_dir"
OBS_DIR="/scratch3/NCEPDEV/marine/Edwin.Rainville/data/sofar_bufr_data"
MODEL_DIR="/scratch3/NCEPDEV/climate/Jessica.Meixner/Data/retrov17_01/"
WW3TOOLS_DIR="/scratch3/NCEPDEV/marine/Edwin.Rainville/ww3-tools-dev"
MODEL_PATTERN="gfs."
MODEL_NAME="GFS_WAVE"
MODEL_FILE_TEMPLATE="gfs.{dt:%Y%m%d}/{dt:%H}/products/wave/gridded/global.0p25/gfs.t{dt:%H}z.global.0p25.f*.grib2"
# ======================================================

# ================= Environment Setup ==================
# Conda 
## Load RDHPCS-conda and activate environment
# ENV_PATH="/scratch3/NCEPDEV/marine/Edwin.Rainville/envs/sofar"
# module load rdhpcs-conda
# eval "$(conda shell.bash hook)"
# conda activate $ENV_PATH

# Prebuilt Modules
module use /scratch4/NCEPDEV/marine/Saeideh.Banihashemi/installs/python-modules/
module load stack-oneapi/2024.2.1
module load ursa-env
module load bufr

# Prioritize the PATH for the URSA environment over the spack stack python
export PATH="/scratch4/NCEPDEV/marine/Saeideh.Banihashemi/installs/ursa-env/bin:$PATH"

# Add the ww3-tools repo root to the python path
PIPELINE="${WW3TOOLS_DIR}/ww3tools/scripts/sofar_matchup_pipeline.py"
export PYTHONPATH="${WW3TOOLS_DIR}${PYTHONPATH:+:${PYTHONPATH}}"

# ======================================================

# Print out the user settings to the log 
echo "==== Submit SOFAR Validation Job ===="
SUBMIT_TIME=$(scontrol show job $SLURM_JOB_ID | grep -o 'SubmitTime=[^ ]*' | cut -d= -f2)
echo "Job Submitted: $SUBMIT_TIME"
echo "Validation Start Date: $START_DATE"
echo "Validation End Date: $END_DATE"
echo "Work Directory: $WORK_DIR"
echo "Observation Directory: $OBS_DIR"
echo "Model Directory: $MODEL_DIR"
echo "WW3tools Directory: $WW3TOOLS_DIR" 
echo "Model directory pattern: $MODEL_PATTERN"
MEM_GB=$((SLURM_MEM_PER_NODE / 1024))
echo "Memory per Node: ${SLURM_MEM_PER_NODE} MB (${MEM_GB} GB)"

# Execute the python matchup script
python -u "$PIPELINE" \
    --start-date "${START_DATE}" \
    --end-date "${END_DATE}" \
    --work-dir "${WORK_DIR}" \
    --obs-dir "$OBS_DIR" \
    --model-dir "${MODEL_DIR}" \
    --model-pattern "${MODEL_PATTERN}" \
    --model-name "${MODEL_NAME}" \
    --model-file-template "${MODEL_FILE_TEMPLATE}"