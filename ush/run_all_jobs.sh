#!/bin/bash

# Script to submit all job_*.sh scripts concurrently in the current directory using sbatch
# Usage: ./run_all_jobs.sh [directory]
# If no directory is provided, uses the current directory

# Set default directory to current if not provided
DIR="${1:-$(pwd)}"

# Maximum number of concurrent jobs (optional, adjust based on Slurm limits)
MAX_CONCURRENT_JOBS=24

# Check if directory exists
if [ ! -d "$DIR" ]; then
  echo "Error: Directory $DIR does not exist."
  exit 1
fi

# Change to the specified directory
cd "$DIR" || { echo "Error: Cannot change to directory $DIR"; exit 1; }

# Log file for tracking submissions
LOGFILE="job_submission_$(date +%Y%m%d_%H%M%S).log"

# Initialize counters
submitted=0
failed=0
queued_jobs=0

echo "Starting job submission from $DIR at $(date)" | tee -a "$LOGFILE"

# Loop through all job_*.sh files
for job_script in job_*.sh; do
  # Check if any job scripts exist
  if [ ! -e "$job_script" ]; then
    echo "No job scripts (job_*.sh) found in $DIR" | tee -a "$LOGFILE"
    exit 0
  fi

  # Check if the file is executable
  if [ ! -x "$job_script" ]; then
    echo "Warning: $job_script is not executable. Making it executable." | tee -a "$LOGFILE"
    chmod +x "$job_script"
  fi

  # Optional: Check number of currently queued/running jobs to avoid overwhelming Slurm
  queued_jobs=$(squeue -u $USER -h | wc -l)
  while [ $queued_jobs -ge $MAX_CONCURRENT_JOBS ]; do
    echo "Too many jobs ($queued_jobs) queued. Waiting 10 seconds..." | tee -a "$LOGFILE"
    sleep 10
    queued_jobs=$(squeue -u $USER -h | wc -l)
  done

  # Submit the job using sbatch
  echo "Submitting $job_script..." | tee -a "$LOGFILE"
  sbatch_output=$(sbatch "$job_script" 2>&1)
  if [ $? -eq 0 ]; then
    echo "Success: $job_script submitted. Output: $sbatch_output" | tee -a "$LOGFILE"
    ((submitted++))
  else
    echo "Error: Failed to submit $job_script. Output: $sbatch_output" | tee -a "$LOGFILE"
    ((failed++))
  fi
done

echo "Submission complete at $(date)" | tee -a "$LOGFILE"
echo "Total jobs submitted: $submitted" | tee -a "$LOGFILE"
echo "Total jobs failed: $failed" | tee -a "$LOGFILE"

if [ $failed -gt 0 ]; then
  echo "Warning: Some jobs failed to submit. Check $LOGFILE for details."
  exit 1
else
  echo "All jobs submitted successfully."
  exit 0
fi
