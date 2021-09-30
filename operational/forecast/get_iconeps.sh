#!/bin/bash

# Script to fetch ICON forecasts
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
PREFIX=${PREFIX:="iconeps"}

# Date of the data to download
DATE=${DATE:=$(date '+%Y%m%d')}
HOUR=${HOUR:="00"}

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

SERVER="https://opendata.dwd.de/weather/nwp/icon-eps/grib"

VARIABLES="u_10m v_10m"

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
[[ ! -d $LOG_DIR ]] && mkdir -p $LOG_DIR
# Init log file
echo 'vim: foldmethod=marker foldlevel=0' > $LOG_FILE

# ======================
# SOURCE FUNCTIONS
source "$(dirname $BASH_SOURCE)/lib/logging.sh"
source "$(dirname $BASH_SOURCE)/lib/util.sh"

# =====================
# FILES USED

WEIGHTS_FILE=${WEIGHTS_FILE:=$(dirname $BASH_SOURCE)/grid/weights_icogl2world_025.nc}
TARGET_FILE=${TARGET_FILE:=$(dirname $BASH_SOURCE)/grid/target_grid_world_025.txt}
GRID_FILE=${GRID_FILE:=$(dirname $BASH_SOURCE)/grid/icon_grid_0024_R02B06_G.nc}

# Get CDO grid description file for global regular grid of ICON.
weight_target_url=https://opendata.dwd.de/weather/lib/cdo/ICON_GLOBAL2WORLD_025_EASY.tar.bz2
grid_url=https://opendata.dwd.de/weather/lib/cdo/icon_grid_0024_R02B06_G.nc.bz2
get_grid_files "$(dirname $BASH_SOURCE)/grid" $weight_target_url $grid_url || exit 1

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

_forecast_hours="$(seq -f '%03g' 0 1 48) $(seq -f '%03g' 51 3 72) $(seq -f '%03g' 78 6 120) $(seq -f '%03g' 132 12 180)"

for var in $VARIABLES; do
	nfiles=0
	for forecast_hour in $_forecast_hours ; do
		[[ $(jobs -p|wc -l) -gt $((${NJOBS}-1)) ]] && wait -n # manage running jobs
		url="${SERVER}/${HOUR}/${var}/icon-eps_global_icosahedral_single-level_${DATE}${HOUR}_${forecast_hour}_${var}.grib2"
		dpath="${_DOWNCAST_CACHE_DIR}/iconeps.${DATE}${HOUR}.${var^^}.f${forecast_hour}"
		download -wget "$url" "$dpath" &&
			post_process -cdo "$dpath" &&
			extract_variables "${dpath}.nc" &
		nfiles=$((nfiles+1))
	done

	wait

	# process all extracted variables (these variables are for each ensemble member 10u_2, 10u_3, 10u_4, 10u_5, etc)
	var_list=$(ls ${_DOWNCAST_CACHE_DIR}/iconeps.${DATE}${HOUR}.${var^^}.f*.nc | awk -F. '{print $(NF-1)}' | sort | uniq)
	downloaded_files=$(ls ${_DOWNCAST_CACHE_DIR}/iconeps.${DATE}${HOUR}.${var^^}.f*.$(echo $var_list | cut -d' ' -f1).nc | wc -l)
	# verify if all files were downloaded, if so concatenate them
	if [[ $nfiles -eq $downloaded_files ]] ; then
		for vv in $var_list ; do
			[[ $(jobs -p|wc -l) -gt $((${NJOBS}-1)) ]] && wait -n # manage running jobs
			var_fmt=$(awk -v vv=$vv 'BEGIN {gsub("10..?","",vv); vv=vv==""?"01":vv; printf "%02d\n", vv}') # convert things like 10u into 01 and 10v_33 into 33
			concatenate_records ${_DOWNCAST_CACHE_DIR}/iconeps.${DATE}${HOUR}.${var^^}.f*.$vv.nc &&
				rename_variables ${_DOWNCAST_CACHE_DIR}/iconeps.${DATE}${HOUR}.${var^^}.$vv.nc ${_RENAME_CONF_FILE} &&
				delete_variables ${_DOWNCAST_CACHE_DIR}/iconeps.${DATE}${HOUR}.${var^^}.$vv.nc "height" &&
				extract_variables ${_DOWNCAST_CACHE_DIR}/iconeps.${DATE}${HOUR}.${var^^}.$vv.nc ${_OUT_DIR}/iconeps.${DATE}${HOUR}.${var_fmt}.nc &
		done
	else
		log_error "$LOG_FILE" "Failed to download/process some files from $var"
	fi
	wait
done

wait

# Cleanup
if [[ $_DOWNCAST_CACHE_DIR != $_OUT_DIR ]] && [[ -d $_DOWNCAST_CACHE_DIR ]] ; then
	rm -r $_DOWNCAST_CACHE_DIR
fi

log "$LOG_FILE" "Finished $(basename $BASH_SOURCE)"

