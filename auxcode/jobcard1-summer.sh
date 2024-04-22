#!/bin/bash

#SBATCH --nodes=1
#SBATCH --partition=orion
#SBATCH --exclusive
#SBATCH --job-name=proc-sat-interpolation
#SBATCH --output=1log-proc_sat_interpolation_%j.out
#SBATCH --time=08:00:00
#SBATCH --account=marine-cpu

# Load necessary modules
module use /work2/noaa/marine/jmeixner/general/modulefiles
module load ww3tools

# Increase stack size
ulimit -s unlimited

# Settings for the job
PYTHON_SCRIPT="ProcSat_interpolation.py"
BASE_MODEL_DATA_DIR="/work/noaa/marine/jmeixner/Data/HR2/summer/"
FILE_FORMAT="grib2"
MODEL="HR2"

# Define start and end datetime for processing
start_datetime="2020060100" 
end_datetime="2020060600"  

# Convert date-times to UNIX timestamps for comparison
start_timestamp=$(date -u -d "${start_datetime:0:8} ${start_datetime:8:2}" +%s)
end_timestamp=$(date -u -d "${end_datetime:0:8} ${end_datetime:8:2}" +%s)

# Iterate over each date and time within the range
current_timestamp=$start_timestamp
while [ $current_timestamp -le $end_timestamp ]; do
    # Format current date and time
    current_datetime=$(date -u -d @$current_timestamp +"%Y%m%d%H")
    current_date=${current_datetime:0:8}
    current_hour=${current_datetime:8:2}


    current_data_dir="${BASE_MODEL_DATA_DIR}gfs.${current_date}/00/products/wave/gridded/" 
    output_data_dir="./"
    # Attempt to find the satellite file within the current data directory
    satellite_file="./Altimeter_CRYOSAT2_HRsummer.nc"

    # Define the pattern for model data files based on the current datetime
    MODEL_DATA_PATTERN='gfswave.t00z.global.0p25.f*.grib2'
    model_data_files=($current_data_dir$MODEL_DATA_PATTERN)

    if [ ${#model_data_files[@]} -gt 0 ]; then
        echo "Processing files for datetime: $current_datetime"
        OUTPUT_FILE="./WW3-Altimeter_summer${MODEL}ainterpolated${current_datetime}.nc"

        # Run the Python script
        python $PYTHON_SCRIPT -t $FILE_FORMAT -d $current_data_dir -p "$MODEL_DATA_PATTERN" -s $satellite_file -o $output_data_dir -f $OUTPUT_FILE -m $MODEL
    else
        echo "No model data files found for datetime: $current_datetime"
    fi

    # Increment current timestamp by 6 hours (86400 seconds)
    current_timestamp=$(($current_timestamp + 86400))
done

