#!/bin/sh
#SBATCH --nodes=1
#SBATCH -q batch
#SBATCH -t 08:00:00
#SBATCH -A marine-cpu
#SBATCH --partition=u1-service
#SBATCH -J get_sofar_hpss
#SBATCH -o get_sofar_hpss.out

# ================= USER'S INPUT =======================
START_DATE="2026-01-01"
END_DATE="2026-01-31"
OBS_DIR="/scratch3/NCEPDEV/marine/Edwin.Rainville/data/sofar_bufr_data"
# ======================================================

# Print the user input to the log
echo "==== Submit SOFAR Data Access Job ===="
SUBMIT_TIME=$(scontrol show job $SLURM_JOB_ID | grep -o 'SubmitTime=[^ ]*' | cut -d= -f2)
echo "Job Submitted: $SUBMIT_TIME"
echo "Data Access Start Date: $START_DATE"
echo "Data Access End Date: $END_DATE"
echo "Observation Directory: $OBS_DIR"

# Move into data directory or give error if it does not exist
cd "$OBS_DIR" || { echo "Error: Cannot cd to $OBS_DIR"; exit 1; }

# Convert dates to seconds for easy looping, adding 16 days to the end date
START_SEC=$(date -d "$START_DATE" +%s)
END_SEC=$(date -d "$END_DATE + 16 days" +%s)

ARCHIVE_BASE="/NCEPPROD/hpssprod/runhistory"
MISSING_FILES=()

# Loop through all days from start date to extended end date
CURRENT_SEC=$START_SEC
while [ $CURRENT_SEC -le $END_SEC ]; do
    YY=$(date -d "@$CURRENT_SEC" +%Y)
    MM=$(date -d "@$CURRENT_SEC" +%m)
    DD=$(date -d "@$CURRENT_SEC" +%d)
    
    TARGET_FILE="sofar_b001_xx121_${YY}${MM}${DD}"
    
    if [ ! -f "$TARGET_FILE" ]; then
        FILE_ATM="${ARCHIVE_BASE}/rh${YY}/${YY}${MM}/${YY}${MM}${DD}/dcom_${YY}${MM}${DD}.tar"
        
        echo "Missing: $TARGET_FILE. Extracting from HPSS at $FILE_ATM ..."
        htar -P -xvf "$FILE_ATM" ./b001
        
        if [ -f "b001/xx121" ]; then
            mv b001/xx121 "$TARGET_FILE"
            echo "Successfully extracted $TARGET_FILE"
        else
            echo "Warning: Target xx121 not found in $FILE_ATM"
            MISSING_FILES+=("$TARGET_FILE")
        fi
        
        # Cleanup the temporary extracted directory
        rm -rf b001
    else
        echo "Found: $TARGET_FILE"
    fi
    
    # Advance exactly 1 day (86400 seconds)
    CURRENT_SEC=$((CURRENT_SEC + 86400))
done

# Final report
if [ ${#MISSING_FILES[@]} -ne 0 ]; then
    echo "=========================================="
    echo "The following files could NOT be extracted:"
    for file in "${MISSING_FILES[@]}"; do
        echo " - $file"
    done
    echo "=========================================="
    exit 1
else
    echo "Success: All requested files are present."
    exit 0
fi
