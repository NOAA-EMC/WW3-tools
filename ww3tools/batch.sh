#!/bin/bash
#SBATCH -q debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:30:00
#SBATCH --partition=hera
#SBATCH --account=marine-cpu
#SBATCH --job-name=process_data
#SBATCH --output=logfile-interpolation.out

# Load modules
module use /scratch1/NCEPDEV/climate/Jessica.Meixner/general/modulefiles
module load module_hera

# Change to the script directory
cd /scratch2/NCEPDEV/marine/Ghazal.Mohammadpour/Tools/WW3-tools/ww3tools/pyresampletest/final/finalist/check/validation/HR1/20200919


# Define paths for the Python script
MODEL_DATA_DIR='/scratch2/NCEPDEV/marine/Jessica.Meixner/Data/HR1/Hurricane/gfs.20200919/00/wave/gridded/'
MODEL_DATA_PATTERN='gfswave.t00z.global.0p25.f*.grib2.nc'
SATELLITE_FILE='./AltimeterAlongTrack_ww3tools_JASON3_2020091901to2020092922.nc'
OUTPUT_FILE='./WW3-Altimeter_interpolated_20200919.nc'

# Run the Python script with paths as arguments
python ProcSat_interpolation.py $MODEL_DATA_DIR $MODEL_DATA_PATTERN $SATELLITE_FILE $OUTPUT_FILE
#python pp.py $MODEL_DATA_DIR $MODEL_DATA_PATTERN $SATELLITE_FILE $OUTPUT_FILE

