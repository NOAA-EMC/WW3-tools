#!/bin/bash
# vim: foldmethod=marker

# Script to perform variable renaming and per variable compression on older files
# -------------------------------------------------------------------------------

# Renaming is automatic, compression (-c) is optional and explicit
#
# There is also the possibility to replace time dimension (-f)
#  for that the argument must take one of the following forms:
#        3hour  6hour 12hour 1day   ...
#
# For files that have 2 different timestep values there is an option (-r) that receives
#  the number of days and a new timestamp.
#  For example if you want : 1h by 1h until the 5th day and 3h by 3h after the 5th day
#         obtaining:  2019-09-10T00:00:00 ...1h... 2019-09-15T00:00:00 ...3h... END
#   the script options would look like this
#        -f 1hour -r 5days 3hour
#      NOTE that if the timesteps need to be fixed with -f option the -r option will not
#            work without them being fixed before (or using -f at the same time)

set -u

# ==============================================================================
# SETTINGS

# Number of parallel jobs to run
NJOBS=${NJOBS:=15}

# ==============================================================================
# COMMAND LINE ARGUMENTS

FILES=()
while [[ $# -gt 0 ]]; do
	arg="$1"

	case $arg in
		-c|--compress)
			COMPRESS=YES
			shift
			;;
		-p|--prefix)
			PREFIX="$2"
			shift # past argument
			shift # past value
			;;
		-f|--fix-record-dim)
			TIME_STEP="$2"
			shift
			shift
			;;
		-r|--rearange-time-steps)
			OFFSET="$2"
			OFFSET_TIME_STEP="$3"
			shift
			shift
			shift
			;;
		*)	  # files
			FILES+=("$1") # save it in an array for later
			shift # past argument
			;;
	esac
done
set -- "${FILES[@]}"

# ======================
# SOURCE FUNCTIONS
source "$(dirname $BASH_SOURCE)/lib/logging.sh"
source "$(dirname $BASH_SOURCE)/lib/util.sh"


# Function to use bellow, transforms record variable into time, and ads a scale to it
## Arguments:
#  $1 -- NETCDF file in which to perform the operations
#  $2 -- a timestep value one of the following:  3hour  6hour 12hour 1day   ...
fix_record_var() {
# {{{
	file=$1
	TIME_STEP=$2
	tmp_filename="/dev/shm/$(basename $file)"
	if [[ $(ncks -mq ${file} | grep -Pzl '(?s)record = UNLIMITED.*\n.*time = 1') ]] ; then
		datetime=$(echo ${file} | sed 's/.*\(20[0-9]\{2\}\)\([0-9]\{2\}\)\([0-9]\{2\}\)\([0-9]\{2\}\).*/\1-\2-\3,\4:00:00/')
		if [[ ${datetime} == ${file} ]] ; then # wasnt able to extract date from filename
			echo "Couldn't extract date from filename ${file}. Impossible to fix record dimension correctly"
			exit 1
		else
			log "$LOG_FILE" "${file}" &&
			{
				ncrename -O -v .time,time_orig -d .time,time_orig ${file} ${tmp_filename}.tmp &&
					ncrename -O -d .record,time ${tmp_filename}.tmp ${tmp_filename} &&
					delete_variables ${tmp_filename} time_orig &&
					cdo -O -z zip -settaxis,${datetime},${TIME_STEP} ${tmp_filename} ${tmp_filename}.tmp &&
					cdo -O -z zip -setreftime,${datetime},hours ${tmp_filename}.tmp ${file} &&
					rm ${tmp_filename}.tmp ${tmp_filename}
			} &>> $LOG_FILE || # append stdout to file, dup stderr to file and terminal
				{ log_error "$LOG_FILE" "Something failed when trying to fix 'time' in $file" && return 1; } # error if one of the previous operation fails
		fi
	fi
# }}}
}


rearange_record_var() {
# {{{
	file=$1
	OFFSET=$2
	OFFSET_TIME_STEP=$3
	tmp_filename1="/dev/shm/$(basename $file).1"
	tmp_filename2="/dev/shm/$(basename $file).2"
	if [[ $(ncks -mq ${file} | grep -Pzl '(?s)record = UNLIMITED.*\n.*time = 1') ]] ; then
		log_error "$LOG_FILE" "Trying to rearange time in $file that hasnt been fixed. Please use -f option." && return 1;
	else
		log "$LOG_FILE" "${file}" && {
			datetime=$(echo ${file} | sed 's/.*\(20[0-9]\{2\}\)\([0-9]\{2\}\)\([0-9]\{2\}\)\([0-9]\{2\}\).*/\1-\2-\3,\4:00:00/')
			first_cuttoff=$(date "+%Y-%m-%dT%H:%M:%S" -d "$( date -d "${datetime/,/T}") +${OFFSET}")
			# in order to remove the final one to avoid duplicates when concatenating
			first_cuttoff_adjusted=$(date "+%Y-%m-%dT%H:%M:%S" -d "$( date -d "${datetime/,/T}") +${OFFSET} -1minute")
			# this is just so we can go until the last timestep (1year should be enough)
			last_timestep=$(date "+%Y-%m-%dT%H:%M:%S" -d "$( date -d "${datetime/,/T}") +1year")
			cdo -O -z zip -seldate,${datetime/,/T},${first_cuttoff_adjusted} $file $tmp_filename1 &&
			cdo -O -z zip -seldate,${first_cuttoff},${last_timestep} $file ${tmp_filename2}.tmp &&
			cdo -O -z zip -settaxis,${first_cuttoff/T/,},${OFFSET_TIME_STEP} ${tmp_filename2}.tmp ${tmp_filename2} &&
			ncrcat -O ${tmp_filename1} ${tmp_filename2} ${file} &&
			rm ${tmp_filename1} ${tmp_filename2} ${tmp_filename2}.tmp
		} &>> $LOG_FILE || # append stdout to file, dup stderr to file and terminal
			{ log_error "$LOG_FILE" "Something failed when trying to rearange 'time' in $file" && return 1; } # error if one of the previous operation fails
	fi
# }}}
}

# ==============================================================================
# COMMAND LINE ARGUMENTS

if [[ "${PREFIX:=NO}" == "NO" ]] ; then # get the prefix from the filename
	echo -e "Can't run without a prefix. Give it with option -p"
	exit 1
fi

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

if [[ -z ${_RENAME_CONF_FILE} ]] ; then
	echo "No rename config file for prefix ${PREFIX}"
	exit 1
fi


# get compress config file
if [[ "${COMPRESS:=NO}" == "YES" ]] && [[ -z ${_COMPRESS_CONF_FILE} ]] ; then
	echo "No compress config file for prefix ${PREFIX}"
	exit 1
fi


for file in $* ; do
	[[ $(jobs -p|wc -l) -gt $((${NJOBS}-1)) ]] && wait -n # manage running jobs

	{

		rename_variables ${file} ${_RENAME_CONF_FILE} &&

		# fix record dimension
		# this process deletes original time dimension, and replaces it with record,
		#  then sets its values correctly base on the filename date
		if [[ ! "${TIME_STEP:=NO}" == "NO" ]] ; then
			fix_record_var ${file} ${TIME_STEP}
		fi

		if [[ ! "${OFFSET:=NO}" == "NO" ]] ; then
			rearange_record_var ${file} ${OFFSET} ${OFFSET_TIME_STEP}
		fi

		if [[ "${COMPRESS:=NO}" == "YES" ]] ; then
			compress_variables ${file} ${file} ${_COMPRESS_CONF_FILE}
		fi

		# remove height variable present in some files
		if [[ $(ncks -mq ${file} | grep 'height = 1') ]] ; then
			delete_variables ${file} height
		fi

	} &
done

wait

