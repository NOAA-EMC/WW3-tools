#!/bin/bash
set -euo pipefail

# =============================================================================
# Auto-detect script directory (important!)
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

# =============================================================================
# Configuration
# =============================================================================
JOB_GLOB="job_*.sh"
BATCH=20
MAX_WAIT_SECONDS=2400   # 40 minutes

# Log file
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${SCRIPT_DIR}/run_all_jobs_${TIMESTAMP}.log"

# =============================================================================
# Logging setup
# =============================================================================
echo "Starting master script at $(date)" | tee "${LOG_FILE}"
echo "Running inside directory: ${SCRIPT_DIR}" | tee -a "${LOG_FILE}"
echo "Job glob: ${JOB_GLOB}" | tee -a "${LOG_FILE}"
echo "Batch size: ${BATCH}" | tee -a "${LOG_FILE}"
echo "Max wait per batch: $((MAX_WAIT_SECONDS/60)) minutes" | tee -a "${LOG_FILE}"
echo "===========================================" | tee -a "${LOG_FILE}"

exec > >(tee -a "${LOG_FILE}") 2>&1

# =============================================================================
# Find job scripts
# =============================================================================
mapfile -t JOB_SCRIPTS < <(ls ${JOB_GLOB} 2>/dev/null | sort -V)

if [[ ${#JOB_SCRIPTS[@]} -eq 0 ]]; then
  echo "ERROR: No job scripts found matching ${JOB_GLOB}" >&2
  exit 1
fi

echo "FOUND ${#JOB_SCRIPTS[@]} job scripts:"
printf '%5s  %s\n' "#" "Filename"
printf '%5s  %s\n' "-----" "------------------------------"
for i in "${!JOB_SCRIPTS[@]}"; do
  printf '%5d  %s\n' $((i+1)) "${JOB_SCRIPTS[$i]}"
done
echo "---------------------------------------"
echo

# =============================================================================
# Submit in batches
# =============================================================================
total=${#JOB_SCRIPTS[@]}
batch_id=1
idx=0

while [[ $idx -lt $total ]]; do
  echo "=================================================="
  echo "Batch ${batch_id}: submitting jobs $((idx+1)) to $((idx+BATCH < total ? idx+BATCH : total))"
  echo "=================================================="

  job_ids=()

  for ((k=0; k<BATCH && idx<total; k++, idx++)); do
    job_script="${JOB_SCRIPTS[$idx]}"

    echo "Submitting: ${job_script}"

    submit_out=$(sbatch --parsable "${job_script}")
    job_id="${submit_out%%;*}"

    if [[ -n "${job_id}" ]]; then
      job_ids+=("${job_id}")
      echo "  -> JobID: ${job_id}"
    else
      echo "WARNING: Failed to parse JobID (${submit_out})"
    fi
  done

  # Wait for batch
  if [[ ${#job_ids[@]} -gt 0 ]]; then
    echo "Waiting for ${#job_ids[@]} jobs to finish..."

    start_time=$(date +%s)
    while true; do
      remaining=$(squeue -h -j "$(IFS=,; echo "${job_ids[*]}")" --states=PENDING,RUNNING 2>/dev/null | wc -l || echo 0)

      if [[ $remaining -eq 0 ]]; then
        echo "Batch ${batch_id} completed."
        break
      fi

      elapsed=$(( $(date +%s) - start_time ))
      if [[ $elapsed -ge $MAX_WAIT_SECONDS ]]; then
        echo "Batch ${batch_id} timeout reached. Continue next batch."
        squeue -h -j "$(IFS=,; echo "${job_ids[*]}")" -o "%i %j %T %R" 2>/dev/null || true
        break
      fi

      sleep 30
    done
  fi

  echo
  batch_id=$((batch_id+1))
done

echo "=================================================="
echo "ALL JOBS SUBMITTED"
echo "Finished at $(date)"
echo "Log: ${LOG_FILE}"
echo "=================================================="
