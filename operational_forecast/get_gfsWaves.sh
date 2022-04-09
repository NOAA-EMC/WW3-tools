#!/bin/bash

# get_gfsWaves.sh
#
# VERSION AND LAST UPDATE:
#  v2.0  10/01/2020
#  v1.1  07/01/2020
#  v1.0  08/01/2017
#
# PURPOSE:
#  Script to download NOAA Global Forecast System (deterministic), Wave 
#   Forecast from WAVEWATCH III. Download from ftp, not nomads.
#  By default the downloads are made into the ./data directory of the repository
#   containing this script, and the logging is made into ./log/<date>/
#  To speed up the download and processing, multiple jobs can be used,
#   taking advantage of multiple cores, if they are available (see and 
#   edit NJOBS below).
#
# USAGE:
#  It can be directly run with ./get_gfsWaves.sh
#   or informing the programs/libraries location (e.g., which wgrib2):
#   WGRIB2=/usr/local/grib2/wgrib2 /home/user/forecast/get_gfsWaves.sh
#  It can also be added to crontab file, so the download is done automaticaly:
#
#   00 06 * * * WGRIB2=/usr/local/grib2/wgrib2 /home/user/forecast/get_gfsWaves.sh
#
#  Customizing the dowload, the variables defined in the next section 
#   #CONTROLABLE ARGUMENTS can be passed to this script, such as
#
#    $ NJOBS=40 DATA_DIR=/my/preferred/directory/ forecast/get_cfs.sh
#
#   or LATMIN, LATMAX, LONMIN, LONMAX
#  Besides the variables defined, it is also possible to define the used
#   programs like wgrib2, cdo, wget, etc. 
#  These programs are defined in the functions file forecast/lib/util.sh
#   however it is still possible to control its values in the same way. 
#   For example, if the script is not finding your installation for 
#   wgrib2 the program path can be given with an uppercase variable with
#   the same name. Ex:
#
#   $ WGRIB2=/usr/local/grib2/wgrib2/wgrib2 forecast/get_cfs.sh
#
#   This can be usefull for running these scripts in a crontab since the
#   path might not contain some programs.
# 
# OUTPUT:
#  One cycle of NOAA's GFS Waves forecast. Netcdf format.
#  A directory is created with format YYYYMMDD00 and one file 
#   ww3.global.0p25.YYYYMMDD00.nc
#
# DEPENDENCIES:
#  wget, perl, CDO, NCO, wgrib2
#
# AUTHOR and DATE:
#  08/01/2017: Ricardo M. Campos, first version based on Ronaldo 
#    Palmeira's codes using perl scripts to select specific variables 
#    from NOAA .idx files.
#  07/01/2020: Ricardo M. Campos, first version for GFS Waves, including 
#    region selection, variables, and post-processing (compression).
#  10/01/2020: Fabio Almeida. New structure and great improvement using
#    shell functions, better coding, new config and lib directory (see util.sh),
#    and option to inform each programs/libraries full path.
#
# PERSON OF CONTACT:
#  Ricardo M Campos: ricardo.campos@noaa.gov
#

set -u

# ==============================================================================
# CONTROLABLE ARGUMENTS
# These arguments can be given through the command line

# Name that is given to the data directory created and the log file
PREFIX=${PREFIX:="gfsWaves"}

# Date of the data to download
DATE=${DATE:=$(date '+%Y%m%d')}
HOUR=${HOUR:="00"}  # cycle

# Number of parallel jobs to run
NJOBS=${NJOBS:=20}

# ===== Specific download variables
# Number of max retries for downloading a file
MAX_RETRIES=${MAX_RETRIES:=60} # retrying ~4:30 hours later
# Minimum size of a file download
MIN_SIZE=${MIN_SIZE:=100000}
# window size to cut the downloaded files
LATMIN=${LATMIN:=-77.5}
LATMAX=${LATMAX:=90.}
LONMIN=${LONMIN:=-102.}
LONMAX=${LONMAX:=30.}

# Data and log directories
DATA_DIR=${DATA_DIR:=$(realpath $(dirname $BASH_SOURCE)/../data)}
[[ ! -d "$(dirname $BASH_SOURCE)/../log" ]] && mkdir -p "$(dirname $BASH_SOURCE)/../log"
LOG_DIR=${LOG_DIR:=$(realpath $(dirname $BASH_SOURCE)/../log/$(date '+%Y%m%d'))}

# ==============================================================================
# DEFINITIONS

SERVER="https://www.ftp.ncep.noaa.gov/data/nccf/com/gfs/prod"

VARIABLES=":UGRD:surface:|:VGRD:surface:|:HTSGW:surface:|:WVHGT:surface:|:SWELL:1 in sequence:|:SWELL:2 in sequence:|:PERPW:surface:|:WVPER:surface:|:SWPER:1 in sequence:|:SWPER:2 in sequence:|:DIRPW:surface:|:WVDIR:surface:|:SWDIR:1 in sequence:|:SWDIR:2 in sequence:"

LOG_FILE="$LOG_DIR/${PREFIX}_${DATE}_${HOUR}.log"

# check for renaming/compress rules files
_RENAME_CONF_FILE=''
_COMPRESS_CONF_FILE=''
_CONF_PATH="$(dirname $BASH_SOURCE)/config/${PREFIX}"
if [[ -d ${_CONF_PATH} ]] ; then
	if [[ -f "$(realpath ${_CONF_PATH}/renames.conf)" ]] ; then
		_RENAME_CONF_FILE="$(realpath ${_CONF_PATH}/renames.conf)"
	fi
	if [[ -f "$(realpath ${_CONF_PATH}/compress.conf)" ]] ; then
		_COMPRESS_CONF_FILE="$(realpath ${_CONF_PATH}/compress.conf)"
	fi
fi

# =======================
# CREATE DIRECTORIES

[[ ! -d $DATA_DIR ]] && mkdir -p $DATA_DIR
[[ ! -d ${DATA_DIR}/${PREFIX}/${DATE}${HOUR} ]] && mkdir -p ${DATA_DIR}/${PREFIX}/${DATE}${HOUR}
[[ ! -d $LOG_DIR ]] && mkdir -p $LOG_DIR
# Init log file
echo 'vim: foldmethod=marker foldlevel=0' > $LOG_FILE

# ======================
# SOURCE FUNCTIONS
source "$(dirname $BASH_SOURCE)/lib/logging.sh"
source "$(dirname $BASH_SOURCE)/lib/util.sh"

# =====================
# OUTPUT DIRECTORY

_OUT_DIR="${DATA_DIR}/${PREFIX}/${DATE}${HOUR}/"
[[ ! -d ${_OUT_DIR} ]] && mkdir -p ${_OUT_DIR}

# =====================
# CACHE
# create a cache dir in memory to make file processing faster
# If there isn't enough space, intermediate files are stored in the output folder
_shm_size=$(df --output=avail /dev/shm | tail -1)
if [[ ${_shm_size} -gt  10000000 ]] ; then # bigger than ~10G
	_DOWNCAST_CACHE_DIR="$(mktemp -d -p /dev/shm/ ${PREFIX}.XXXXX || echo ${_OUT_DIR})"
else
	_DOWNCAST_CACHE_DIR="${_OUT_DIR}"
fi

# ============================
# MAIN

_forecast_hours="$(seq -f '%03g' 0 1 120) $(seq -f "%03g" 123 3 384)"

nfiles=0
for forecast_hour in $_forecast_hours ; do
	[[ $(jobs -p|wc -l) -gt $((${NJOBS}-1)) ]] && wait -n # manage running jobs
	url="${SERVER}/gfs.${DATE}/${HOUR}/wave/gridded/gfswave.t${HOUR}z.global.0p25.f${forecast_hour}.grib2"
	dpath="${_DOWNCAST_CACHE_DIR}/${PREFIX}.global.0p25.${DATE}${HOUR}.f${forecast_hour}"
	download -get-perl "$url" "$dpath" "$VARIABLES" &&
		post_process -wgrib "$dpath" &
	nfiles=$((nfiles+1))
done

wait
file_list=${_DOWNCAST_CACHE_DIR}/${PREFIX}.global.0p25.${DATE}${HOUR}.f*.nc
downloaded_files=$(echo $file_list | tr ' ' '\n' | wc -l)
filename="${PREFIX}.global.0p25.${DATE}${HOUR}.nc"
# verify if all files were downloaded, if so concatenate them
if [[ $nfiles -eq $downloaded_files ]] ; then
	concatenate_records ${_DOWNCAST_CACHE_DIR}/${filename} -- $file_list &&
		rename_variables ${_DOWNCAST_CACHE_DIR}/${filename} ${_RENAME_CONF_FILE} &&
		compress_variables ${_DOWNCAST_CACHE_DIR}/${filename} ${_OUT_DIR}/${filename} ${_COMPRESS_CONF_FILE} &
else
	log_error "$LOG_FILE" "Failed to download/process some files"
fi

wait

# Cleanup
if [[ $_DOWNCAST_CACHE_DIR != $_OUT_DIR ]] && [[ -d $_DOWNCAST_CACHE_DIR ]] ; then
	rm -r $_DOWNCAST_CACHE_DIR
fi

log "$LOG_FILE" "Finished $(basename $BASH_SOURCE)"

