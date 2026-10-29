#!/bin/bash

# ================= USER'S INPUT =======================
START_DATE="2026-01-01"
END_DATE="2026-01-31"
OBS_DIR="/scratch3/NCEPDEV/marine/Edwin.Rainville/data/sofar_bufr_data"
# ======================================================

echo "==== Submitting Parallel SOFAR HPSS Jobs ===="
echo "Data Access Start Date: $START_DATE"
echo "Data Access End Date: $END_DATE"
echo "Observation Directory: $OBS_DIR"

# Move into data directory or give error if it does not exist
mkdir -p "$OBS_DIR"
cd "$OBS_DIR" || { echo "Error: Cannot cd to $OBS_DIR"; exit 1; }

# Convert dates to seconds for easy looping, adding 16 days to the end date
START_SEC=$(date -d "$START_DATE" +%s)
END_SEC=$(date -d "$END_DATE + 16 days" +%s)

ARCHIVE_BASE="/NCEPPROD/hpssprod/runhistory"

# Loop through all days from start date to extended end date
CURRENT_SEC=$START_SEC
while [ $CURRENT_SEC -le $END_SEC ]; do
    YY=$(date -d "@$CURRENT_SEC" +%Y)
    MM=$(date -d "@$CURRENT_SEC" +%m)
    DD=$(date -d "@$CURRENT_SEC" +%d)
    
    TARGET_FILE="sofar_b001_xx121_${YY}${MM}${DD}"
    
    if [ ! -f "$TARGET_FILE" ]; then
        echo "Missing: $TARGET_FILE. Submitting extraction job..."
        
        # Submit a unique job for this specific day using a Here-Doc
        sbatch << EOF
#!/bin/sh
#SBATCH --nodes=1
#SBATCH -q batch
#SBATCH -t 08:00:00
#SBATCH -A marine-cpu
#SBATCH --partition=u1-service
#SBATCH -J sofar_${YY}${MM}${DD}
#SBATCH -o ${OBS_DIR}/sofar_hpss_${YY}${MM}${DD}_%j.out

cd "${OBS_DIR}"

# Create a unique temporary directory to prevent race conditions with parallel jobs
TMP_DIR="tmp_extract_${YY}${MM}${DD}_\$SLURM_JOB_ID"
mkdir -p "\$TMP_DIR"
cd "\$TMP_DIR"

FILE_ATM="${ARCHIVE_BASE}/rh${YY}/${YY}${MM}/${YY}${MM}${DD}/dcom_${YY}${MM}${DD}.tar"
echo "Extracting from HPSS at \$FILE_ATM ..."

# Extract into the unique temporary directory
htar -P -xvf "\$FILE_ATM" ./b001

if [ -f "b001/xx121" ]; then
    # Move the extracted file back up to the main observation directory
    mv b001/xx121 "${OBS_DIR}/${TARGET_FILE}"
    echo "Successfully extracted ${TARGET_FILE}"
else
    echo "Warning: Target xx121 not found in \$FILE_ATM"
fi

# Clean up temporary directory
cd "${OBS_DIR}"
rm -rf "\$TMP_DIR"
EOF
    else
        echo "Found: $TARGET_FILE (Skipping)"
    fi
    
    # Advance exactly 1 day (86400 seconds)
    CURRENT_SEC=$((CURRENT_SEC + 86400))
done

echo "==== Job submission complete. ===="
echo "Monitor your jobs using 'squeue -u \$USER'"