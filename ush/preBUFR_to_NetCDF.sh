#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH -q batch
#SBATCH -t 08:00:00
#SBATCH -A marine-cpu
#SBATCH -J jason3_pb2nc
#SBATCH -o run_JASON3_pb2nc.o%j
#SBATCH -e run_JASON3_pb2nc.e%j

## This script used to convert preBUFR JASON3 satellite data to NetCDF
## It will first discovering all available preBUFR files
## Then check the netCDF DIR, converting only newly added ones (skipping already converted dates)

set -euo pipefail

# ================= USER'S INPUT =======================
MACHINE="ursa"
WORKDIR="/scratch3/NCEPDEV/marine/Ming.Chen/ursa/wind_eval"            # working directory contains WW3-tools

IN_DIR="/scratch3/NCEPDEV/marine/Saeideh.Banihashemi/RWPS/WND/JASON3"  # preBUFR files directory

# ======================================================

OUT_DIR="${WORKDIR}/processsatdata/pb2nc_out"                          #path to pb2nc_nc
CFG="${WORKDIR}/WW3-tools/parm/jason3_config.conf"                     #path to jason3_config.conf

# load necessary modules including METplus (will add hercules and orion soon)
if [[ "$MACHINE" == "ursa" ]]; then
    module use /scratch4/NCEPDEV/marine/Saeideh.Banihashemi/installs/python-modules/
    module load ursa-env
    module load met/12.0.1
elif [[ "$MACHINE" == "hercules" ]]; then
    echo "ERROR: No module configuration defined for machine $MACHINE"
    exit 1
else
    echo "ERROR: No module configuration defined for machine $MACHINE"
    exit 1
fi

mkdir -p "$OUT_DIR"

log() {
  echo "[$(date -u +'%Y-%m-%dT%H:%M:%SZ')] $*"
}

log "===== pb2nc incremental run started ====="
log "Host: $(hostname)"
log "PWD:  $(pwd)"
log "IN_DIR:  $IN_DIR"
log "OUT_DIR: $OUT_DIR"
log "CFG:     $CFG"

# check pb2nc configuration file
if [[ ! -r "$CFG" ]]; then
  log "ERROR: config file not readable: $CFG"
  exit 2
fi

# check all preBUFR files
mapfile -t inputs < <(
  find "$IN_DIR" -maxdepth 1 -type f \
    -name "jason3_b031_xx124_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]" | sort
)

if (( ${#inputs[@]} == 0 )); then
  log "No input files found. Exiting."
  exit 0
fi

log "Found ${#inputs[@]} total preBUFR files."

# Selet only unconverted files
to_do=()
skipped=0

for in_file in "${inputs[@]}"; do
  base=$(basename "$in_file")
  ymd="${base##*_}"
  out_file="${OUT_DIR}/jason3_b031_xx124_${ymd}.nc"

  if [[ -s "$out_file" ]]; then
    ((++skipped))
  else
    to_do+=("$in_file")
  fi
done

log "Already converted (skipped): $skipped"
log "New files to convert:        ${#to_do[@]}"

if (( ${#to_do[@]} == 0 )); then
  log "Nothing new to convert. Done."
  exit 0
fi

log "Will convert these dates:"
for f in "${to_do[@]}"; do
  b=$(basename "$f"); log "  ${b##*_}"
done

log "Starting pb2nc conversions..."

count=0

for in_file in "${to_do[@]}"; do
  base=$(basename "$in_file")
  ymd="${base##*_}"
  out_file="${OUT_DIR}/jason3_b031_xx124_${ymd}.nc"

  log "CONVERT $ymd  in=$in_file  out=$out_file"

  pb2nc \
    "$in_file"  \
    "$out_file" \
    "$CFG"      \
    -v 2

  ((++count))

  # Throttle to at most SLURM_NTASKS concurrent tasks
  if (( count % SLURM_NTASKS == 0 )); then
    wait
  fi
done

  # Wait for any remaining tasks
wait

log "All pb2nc conversions finished."
