#!/bin/bash

# Script to fetch NCEP Ensemble forecasts GDAS
# --------------------------------------------
# Gets forecast of the previous day.
#
# By default the downloads are made into the ./data directory of the repository
#  containing this script, and the logging is made into ./log/<date>/
#
# The variables defined in the next section #CONTROLABLE ARGUMENTS can be passed
#  to this script by defining them when calling the script, ex:
#
#      $ NJOBS=40 DATA_DIR=/my/preferred/directory/ forecast/get_cfs.sh
#
# Besides the variables defined in this section it is also possible to define
#  the used programs like wgrib2, cdo, wget, etc. , these programs are defined
#  in the functions file forecast/lib/util.sh however it is still
#  possible to control its values in the same way. For example if the script is
#  not finding your installation for wgrib2 the program path can be given with
#  an uppercase variable with the same name. Ex:
#
#      $ WGRIB2=/usr/local/grib2/wgrib2/wgrib2 forecast/get_cfs.sh
#
#  This can be usefull for running these scripts in a crontab since the path might
#   not contain some programs.

set -u

# ==============================================================================
# CONTROLABLE ARGUMENTS
# These arguments can be given through the command line

# Name that is given to the data directory created and the log file
PREFIX=${PREFIX:="gdas"}

# Date of the download
DATE=${DATE:=$(date '+%Y%m%d')}
# Number of days to download (before date defined above)
NDAYS=${NDAYS:=4}

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

VARIABLES=":UGRD:850 mb:|:UGRD:10 m above ground:|:VGRD:850 mb:|:VGRD:10 m above ground:|:GUST:surface:|:HGT:850 mb:|:ICEC:surface:|:RH:2 m above ground:|:TMP:2 m above ground:|:PRMSL:mean sea level:"

LOG_FILE="$LOG_DIR/${PREFIX}_${DATE}.log"

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
[[ ! -d $LOG_DIR ]] && mkdir -p $LOG_DIR
# Init log file
echo 'vim: foldmethod=marker foldlevel=0' > $LOG_FILE

# ======================
# SOURCE FUNCTIONS
source "$(dirname $BASH_SOURCE)/lib/logging.sh"
source "$(dirname $BASH_SOURCE)/lib/util.sh"

# =====================
# OUTPUT DIRECTORY

_OUT_DIR="${DATA_DIR}/${PREFIX}/${DATE}00/"
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

_cycles=$(seq -f '%02g' 0 6 18)

nfiles=0
for day in $(seq $NDAYS -1 1) ; do
	nth_date=$(date -u --date="${DATE} -$day day" '+%Y%m%d')

	for cycle in $_cycles ; do

		if [[ ${cycle} -eq "18" ]] && [[ ${day} -eq "001" ]] ; then
			_forecast_hours="$(seq -f '%03g' 0 1 6)"
		else
			_forecast_hours="$(seq -f '%03g' 0 1 5)"
		fi

		for forecast_hour in $_forecast_hours; do
			[[ $(jobs -p|wc -l) -gt $((${NJOBS}-1)) ]] && wait -n # manage running jobs
			url="${SERVER}/gdas.${nth_date}/${cycle}/atmos/gdas.t${cycle}z.pgrb2.0p25.f${forecast_hour}"
			dpath="${_DOWNCAST_CACHE_DIR}/gdas.${nth_date}.t${cycle}z.pgrb2a.0p25.f${forecast_hour}"
			download -get-perl "$url" "$dpath" "$VARIABLES" &&
				post_process -wgrib "$dpath" &
			nfiles=$((nfiles+1))
		done
	done
done

wait
file_list=${_DOWNCAST_CACHE_DIR}/gdas.*.f*.nc
downloaded_files=$(echo $file_list | tr ' ' '\n' | wc -l)
filename_nowcast="gdas.${DATE}.nowcast.pgrb2a.0p25.nc"
filename="gdas.${DATE}.last${NDAYS}days.pgrb2a.0p25.nc"
# verify if all files were downloaded, if so concatenate them
if [[ $nfiles -eq $downloaded_files ]] ; then
	# save nowcast file
	mv -bf "${_DOWNCAST_CACHE_DIR}/gdas.$(date --date="${DATE} -1 day" '+%Y%m%d').t18z.pgrb2a.0p25.f006.nc" ${_DOWNCAST_CACHE_DIR}/${filename_nowcast}
	# rename and compress nowcast
	rename_variables ${_DOWNCAST_CACHE_DIR}/${filename_nowcast} ${_RENAME_CONF_FILE} &&
		compress_variables ${_DOWNCAST_CACHE_DIR}/${filename_nowcast} ${_OUT_DIR}/${filename_nowcast} ${_COMPRESS_CONF_FILE} &
	# concatenate the rest and rename the output file
	concatenate_records ${_DOWNCAST_CACHE_DIR}/${filename} -- $file_list &&
		rename_variables ${_DOWNCAST_CACHE_DIR}/${filename} ${_RENAME_CONF_FILE} &&
		compress_variables ${_DOWNCAST_CACHE_DIR}/${filename} ${_OUT_DIR}/${filename} ${_COMPRESS_CONF_FILE} &
else
	log_error "$LOG_FILE" "Failed to download/process some files from ${DATE}"
fi

wait

# Cleanup
if [[ $_DOWNCAST_CACHE_DIR != $_OUT_DIR ]] && [[ -d $_DOWNCAST_CACHE_DIR ]] ; then
	rm -r $_DOWNCAST_CACHE_DIR
fi

log "$LOG_FILE" "Finished $(basename $BASH_SOURCE)"

