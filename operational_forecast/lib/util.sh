# vim: foldmethod=indent foldnestmax=1

# This file contains functions meant to be used by other scripts
# NOTE:
# The functions defined here expect a few variables to be set:
#   LOG_FILE -- This variable must contain the name of the main log file
#                if nothing is given this will default to a temporary file in /tmp/downcast.XXXXX.log
#
# 10/01/2020: Fabio Almeida. 

# Here is a list of the dependencies of these functions, if they don't exist in the PATH then
#  a variable with the program name in uppercase can be defined with the command to execute for each.
#  For example, cdo can be defined as CDO=/usr/local/bin/cdo in case the script can't find cdo.
#
# ncatted
# ncks
# ncrcat
# ncrename
# ncwa
# wget
# wgrib2
# cdo

source "$(dirname $BASH_SOURCE)/logging.sh"

set -u
set -o pipefail

# Number of max retries for downloading a file
MAX_RETRIES=${MAX_RETRIES:=60} # retrying ~4:30 hours later

# Minimum size of a file download
MIN_SIZE=${MIN_SIZE:=100000}

# window size to cut the downloaded files
LATMIN=${LATMIN:=-77.5}
LATMAX=${LATMAX:=90.}
LONMIN=${LONMIN:=-102.}
LONMAX=${LONMAX:=30.}

# set LOG_FILE to some tmpfile in case nothing is given
LOG_FILE=${LOG_FILE:="$(mktemp /tmp/downcast_log.XXXXX)"}

# verify existance of necessary scrips
[[ -x ${GET_GRIB:=$(realpath $(dirname $BASH_SOURCE)/get_grib.pl)} ]] || { log_error "$LOG_FILE" "get_grib.pl does not exist. A custom path can be given to the program with the GET_GRIB variable." ; exit 1 ;}
[[ -x ${GET_INV:=$(realpath $(dirname $BASH_SOURCE)/get_inv.pl)} ]]   || { log_error "$LOG_FILE" "get_inv.pl does not exist. A custom path can be given to the program with the GET_INV variable." ; exit 1 ;}

# Verify existance of all the necessary programs
[[ -x ${NCATTED:=$(which ncatted)} ]]   || { log_error "$LOG_FILE" "ncatted does not exist. A custom path can be given to the program with the NCATTED variable." ; exit 1 ;}
[[ -x ${NCKS:=$(which ncks)} ]]         || { log_error "$LOG_FILE" "ncks does not exist. A custom path can be given to the program with the NCKS variable." ; exit 1 ;}
[[ -x ${NCRCAT:=$(which ncrcat)} ]]     || { log_error "$LOG_FILE" "ncrcat does not exist. A custom path can be given to the program with the NCRCAT variable." ; exit 1 ;}
[[ -x ${NCECAT:=$(which ncecat)} ]]     || { log_error "$LOG_FILE" "ncecat does not exist. A custom path can be given to the program with the NCECAT variable." ; exit 1 ;}
[[ -x ${NCRENAME:=$(which ncrename)} ]] || { log_error "$LOG_FILE" "ncrename does not exist. A custom path can be given to the program with the NCRENAME variable." ; exit 1 ;}
[[ -x ${NCWA:=$(which ncwa)} ]]         || { log_error "$LOG_FILE" "ncwa does not exist. A custom path can be given to the program with the NCWA variable." ; exit 1 ;}
[[ -x ${WGET:=$(which wget)} ]]         || { log_error "$LOG_FILE" "wget does not exist. A custom path can be given to the program with the WGET variable." ; exit 1 ;}
[[ -x ${WGRIB2:=$(which wgrib2)} ]]     || { log_error "$LOG_FILE" "wgrib2 does not exist. A custom path can be given to the program with the WGRIB2 variable." ; exit 1 ;}
# NOTE: CDO must be built with netcdf4 and eccodes support (can be configured like so: ./configure --with-netcdf=/usr/ --with-szlib=/usr/ --with-hdf5=/usr --with-eccodes=/usr/) if all libraries are in /usr
[[ -x ${CDO:=$(which cdo)} ]]           || { log_error "$LOG_FILE" "cdo does not exist. A custom path can be given to the program with the CDO variable." ; exit 1 ;}

# This function downloads a file from a url with 2 download modes
# ---------------------------------------------------------------------
# The download mode defines the way the file is download,
#   -wget     simply runs wget for a given url adding .bz2 to it
#   -get-perl using get scripts in lib/, with wich variables names can be given separated by | to select which variables to download
## Arguments:
#   $1 -- Download mode -wget or -get-perl
#   $2 -- Base url (this url must not contain .idx extension for the inventory, or .bz2 extension for compressed files)
#   $3 -- Filename without extension
#   $4 -- Variable Names (only when $1 == -get-perl)
download() {
	[[ $# -lt 3 ]] && { log_error "$LOG_FILE" "function must be given at least 3 arguments." ; exit 1 ;}
	[[ "$1" == "-get-perl" ]] && [[ $# == 3 ]] && { log_error "$LOG_FILE" "function. Variables must be given as 4th argument when using -get-perl option." ; exit 1 ;}
	local download_mode=$1 url=$2 file=$3 # rename arguments

	# log into custom log file since this code might run concurrently
	local custom_log_file="${LOG_FILE}.$(basename $file).log"

	# Try to download until file exists and has acceptable size
	local file_size=0 tries=1
	while [[ $file_size -lt $MIN_SIZE ]] && [[ $tries -le ${MAX_RETRIES} ]] ; do
		[[ $file_size -gt 0 ]] && [[ $file_size -lt $MIN_SIZE ]] && log_warn "$custom_log_file" "File $file did not reach minimum allowed size. [$file_size/$MIN_SIZE]"
		[[ $tries -gt 5 ]] && sleep 300
		log "$custom_log_file" "$url -- Attempt #${tries}"

		echo '{{{' >> $custom_log_file # open folds in logfile

		# select download mode
		case "$download_mode" in
			-wget) # use wget
				output_file="${file}.grb2.bz2"
				$WGET --no-verbose "${url}.bz2" -O "${output_file}" &>>$custom_log_file
				;;
			-get-perl) # use perl scripts
				output_file="${file}.grb2"
				[[ $# -ne 4 ]] && { log_error "$LOG_FILE" "function must be given 4 arguments with option -get-perl." ; exit 1 ;}
				$GET_INV "${url}.idx" 2>>$custom_log_file | egrep "$4" | $GET_GRIB $url "${output_file}" &>>$custom_log_file
				;;
			*) # unknown option, throw error and stop script
				echo '}}}' >> $custom_log_file # close folds in logfile
				log_error "$custom_log_file" "Unknown download mode given to function ($download_mode). Must be -wget or -get-perl"
				cat $custom_log_file >> $LOG_FILE && rm $custom_log_file
				exit 1
				;;
		esac

		echo '}}}' >> $custom_log_file # close folds in logfile

		# calculate file size
		if [[ -f ${output_file} ]] ; then
			file_size=$(du -sb "${output_file}" | awk '{print $1}')
		else
			file_size=0
			log_warn "$custom_log_file" "Failed to download ${output_file}"
		fi
		tries=$((${tries}+1))
	done

	# Exit in case MAX_RETRIES reached
	if [[ $tries -ge $MAX_RETRIES ]] ; then
		log_error "$custom_log_file" "Reached maximum number of retries for $url"
		cat $custom_log_file >> $LOG_FILE && rm $custom_log_file
		exit 1
	fi

	cat $custom_log_file >> $LOG_FILE && rm $custom_log_file
}

# This function performs multiple post process operations on a given filename
# ---------------------------------------------------------------------
# The process mode defines the way the file is download,
#   -wgrib    simply runs $WGRIB2 to convert .grib2 to .nc file
#   -cdo      uses $CDO and the variables TARGET_FILE and WEIGHTS_FILE to convert a grib2 file into NETCDF
#   -skip     skip conversion to NETCDF. Useful if you have a .nc and just want to perform the remaining operations
# This functions performs the following operations in sequence:
#   - decompress .bz2 if necessary
#   - converts to .nc
#   - moves coordinates (longitude to -180,180 and latitude to -90,90)
#   - renames lat to latitude if necessary (same for longitude)
#   - cuts latitude and longitude innto the values defined by the variables LATMIN, LATMAX, LONMIN, LONMAX
#   - sets the attribute _FillValue to NaN
# ---------------------------------------------------------------------
## Arguments:
#   $1 -- Process mode -wgrib, -cdo or -skip
#   $2 -- Filename without extension
post_process() {
	[[ $# -ne 2 ]] && { log_error "$LOG_FILE" "function must be given 2 arguments." ; exit 1 ;}

	local process_mode=$1 file=$2 # rename arguments
	local custom_log_file="${LOG_FILE}.$(basename $file).log"
	local nc_file="${file}.nc"

	log "$custom_log_file" "$file"
	echo '{{{' >> $custom_log_file # add folds to logfile

	if ! { # do post processing in a single pipeline and send all output to $custom_log_file
		_decompress "$file" "$custom_log_file" &&
			{ if [[ "$process_mode" != "-skip" ]] ; then _convert_to_nc "$process_mode" "$file" "$custom_log_file" ; fi ; } &&
			$CDO -s sellonlatbox,-180,180,-90,90 "$nc_file" "${nc_file}.tmp" &&
			$NCRENAME -v .lat,latitude -v .lon,longitude -d .lat,latitude -d .lon,longitude "${nc_file}.tmp" &&
			$NCKS -O -4 -L 1 -d latitude,${LATMIN},${LATMAX} "${nc_file}.tmp" "${nc_file}.tmp" &&
			$NCKS -O -4 -L 1 -d longitude,${LONMIN},${LONMAX} "${nc_file}.tmp" "${nc_file}" &&
			$NCATTED -a _FillValue,,o,f,NaN "${nc_file}" ;
	} 2> >(tee -a $custom_log_file) 1>>$custom_log_file ; then # append stdout to file, dup stderr to file and terminal
		echo '}}}' >> $custom_log_file # close folds in logfile
		log_error "$custom_log_file" "Failed processing $nc_file"
		cat $custom_log_file >> $LOG_FILE && rm $custom_log_file
		exit 1
	fi

	# Cleanup
	echo '}}}' >> $custom_log_file # close folds in logfile
	rm ${nc_file}*.tmp ${file}.grb2 &>>$custom_log_file
	cat $custom_log_file >> $LOG_FILE && rm $custom_log_file
}

# Download and generate grid files to use when converting .grib2 to .nc with $CDO
# ---------------
# *NOTE* This function expects TARGET_FILE, WEIGHTS_FILE, and GRID_FILE to be defined variables with file paths
## Arguments
#  $1 -- Output directory
#  $2 -- ICON_GLOBAL tar.bz2 url
#  $3 -- grid file url
get_grid_files() {
	[[ $# -ne 3 ]] && { log_error "$LOG_FILE" "function must be given 3 arguments." ; exit 1 ;}

	[[ ! -d $1 ]] && mkdir -p "$1"

	if [[ ! -f ${TARGET_FILE} ]] ; then
		log "$LOG_FILE" "Weights and Target files doesn't exist. Downloading and extracting them into ${1}..."
		{
			$WGET --no-verbose $2 -O ${1}/ICON_GLOBAL2WORLD.tar.bz2 &&
			tar xjf ${1}/ICON_GLOBAL2WORLD.tar.bz2 --directory ${1}/ &&
			mv ${1}/ICON_GLOBAL2*/*.txt ${TARGET_FILE} &&
			mv ${1}/ICON_GLOBAL2*/*.nc ${WEIGHTS_FILE} &&
			rm -r ${1}/ICON_GLOBAL2* ;
		} 2> >(tee -a $LOG_FILE) 1>>$LOG_FILE || # append stdout to file, dup stderr to file and terminal
			{ log_error "$LOG_FILE" "Failed to download grid files. Cannot continue" && return 1; } # error if one of the previous operation fails
	fi
	if [[ ! -f ${GRID_FILE} ]] ; then
		log "$LOG_FILE" "Grid file doesn't exist. Downloading and extracting grid file ${1}..."
		{
			$WGET --no-verbose $3 -O ${GRID_FILE}.bz2 &&
			$BUNZIP ${GRID_FILE}.bz2 ;
		} 2> >(tee -a $LOG_FILE) 1>>$LOG_FILE || # append stdout to file, dup stderr to file and terminal
			{ log_error "$LOG_FILE" "Failed to download grid files. Cannot continue" && return 1; }
	fi

	log "$LOG_FILE" "Generating file with weights for interpolation in ${1}"
	$CDO gennn,${TARGET_FILE} ${GRID_FILE} ${WEIGHTS_FILE} &>>$LOG_FILE
}

# Concatenate multiple forecast files into a single NETCDF.
# ---------------------------------------------------------
# The function deletes all the files given.
## Arguments:
#  $@ -- Multiple NETCDF files to be concatenated (filename must end in .f###.nc)
## OR
#  $1     -- Name of the output file
#  $2     -- '--' literally two dashes to seperate output file from inputs
#  ${@:2} -- Remaining arguments are input files
concatenate_records() {
	shopt -s extglob # enable the parameter expansion expression bellow

	if [[ $2 == '--' ]] ; then
		output_file="${1}"
		shift ; shift # discard first two arguments
	else
		output_file="${1/.f*([^.])./.}"
	fi

	log "$LOG_FILE" "${1%.f*}.f###.nc -> $output_file"

	if ! {
		# when generating the name remove only the forecast time part "f????."
		$NCRCAT $@ -O $output_file &&
			rm $@
	} 2> >(tee -a $LOG_FILE) 1>>$LOG_FILE ; then # append stdout to file, dup stderr to file and terminal
		log_error "$LOG_FILE" "Failed concatenating into ${output_file}" && return 1
	fi
}

# List the variables with more than one dimension for the given netcdf files
list_variables() {
	local f
	for f in "$@" ; do
		{ $NCKS -mq "$f" | grep 'variables:' -A100 | grep '(.*,.*)' | awk '{gsub("\\(.*", "", $2);gsub("\\\\","",$2);print $2}' ; } 2>>$LOG_FILE
	done
}

# Extract each data variable NETCDF into a separate file, if a second argument is given all variables will be appended to this file.
# ----------------------------------------------------------------------------------------------------------------------------------
# In case the second argument is not given the output files will be named as such (assuming input is "filename.nc") ==> "filename.varname.nc"
## Arguments:
#  $1 -- Original NETCDF file from which to extract variables (this file will be removed)
#  $2 -- Destination NETCDF file where to append the variable
extract_variables() {
	[[ $# -lt 1 ]] && { log_error "$LOG_FILE" "function must be given at least 1 argument." ; exit 1 ;}

	local filename="$1"
	local v_list=$(list_variables "$filename")

	log "$LOG_FILE" "from $filename"

	if ! {
		local vname
		for vname in $v_list ; do
			# define output filename
			[[ $# -lt 2 ]] && output_file="${filename%.nc}.$vname.nc" || output_file="$2"
			# if a second argument is given use it as output, otherwise generate the filename
			$NCKS -A -v $vname "$filename" "$output_file"
		done &&
			rm $filename
	} 2> >(tee -a $LOG_FILE) 1>>$LOG_FILE ; then # append stdout to file, dup stderr to file and terminal
		log_error "$LOG_FILE" "Failed extracting variables from $filename" && return 1;
	fi
}

# Renames variables in a NETCDF file given a variable mapping
# -----------------------------------------------------------
# It expects 2 arguments, the file in which to rename the variables and another
#  that contains the variable mappings with a format as follows:
#                # lines with hashtag or that start with space are ignored
#                # variable names cannot contain spaces
#                old_name1 => new_name1
#                old_name2 => new_name2
# See config_renaming/README.md for more information
## Arguments:
#  $1 -- NETCDF file in which to perform renaming
#  $2 -- Config file with variable mappings
rename_variables() {
	if [[ $# == 1 ]] ; then
		log_warn "$LOG_FILE" "No config file exists, renaming not performed."
		return
	else [[ $# != 2 ]] && { log_error "$LOG_FILE" "function must be given 2 values." ; exit 1 ;}
	fi

	# rename inputs
	local filename="$1"
	local conf_file="$2"

	[[ ! -f "$conf_file" ]] && { log_warn "$LOG_FILE" "$conf_file does not exist, renaming not performed"; return; }

	local variables=( $(list_variables $filename) )
	local rename_args=""

	# read config file line by line
	local orig sign target
	while read orig sign target ; do
		# checking if name is contained in variable array with delimiting spaces can be considered safe since
		# variable do not contain spaces and if for some reason they did, the script would have broke long before :)
		if [[ $sign == "=>" ]] && [[ " ${variables[@]} " =~ " $orig " ]] ; then
			rename_args="${rename_args} -v .$orig,$target"
		fi
	done <<< $(grep -v "^\ *#\|^\ *$" "$conf_file")

	# if there is something to rename, rename it
	if [[ ! -z ${rename_args} ]] ; then
		log "$LOG_FILE" "$filename | args:$rename_args"
		$NCRENAME $rename_args "${filename}" &>> $LOG_FILE ||
			{ log_error "$LOG_FILE" "Failed to rename variables in $filename" && return 1; }
	else
		log "$LOG_FILE" "No variables in $conf_file match. Renaming not performed."
		return
	fi
}

# Delete a set of dimensions from a file (degenerate dimensions, length 1)
# ------------------------------------------------------------------------
# This function will use the element zero of a dimension in the data variables, in order to keep
## Arguments:
#  $1     -- NETCDF file in which to perform deletion
#  ${@:2} -- Remaining arguments are variables to delete
delete_variables() {
	local filename="$1"
	shift
	local target_vars="${@}"

	log "$LOG_FILE" "'${@}' in $filename"
	local var
	for var in $target_vars ; do
		{
			$NCKS -4 -C -O -x -v $var "${filename}" "${filename}.tmp" &&
			$NCWA -4 -O -d ${var},0 -a $var "${filename}.tmp" "${filename}" &&
			rm "${filename}.tmp" ;
		} 2> >(tee -a $LOG_FILE) 1>>$LOG_FILE || # append stdout to file, dup stderr to file and terminal
			{ log_error "$LOG_FILE" "Failed to delete $var in $filename" && return 1; } # error if one of the previous operation fails
	done
}


# Compress variables in a NETCDF file given a variable mapping
# -----------------------------------------------------------
# It expects 2 arguments, the file in which to compress the variables and another
#  that contains the variable mappings with a format as follows:
#                # lines with hashtag or that start with space are ignored
#                # variable names cannot contain spaces
#                var_name1 => #_digits_to_keep
#                var_name2 => #_digits_to_keep
# See config_compress/README.md for more information
## Arguments:
#  $1 -- NETCDF file in which to perform compression
#  $2 -- Output NETCDF file
#  $3 -- Config file with variable mappings
compress_variables() {
	if [[ $# == 1 ]] ; then
		log_warn "$LOG_FILE" "No config file exists, compression not performed."
		return
	else [[ $# != 3 ]] && { log_error "$LOG_FILE" "function must be given 3 values." ; exit 1 ;}
	fi

	# rename inputs
	local filename="$1"
	local output_file="$2"
	local conf_file="$3"

	[[ ! -f "$conf_file" ]] && { log_warn "$LOG_FILE" "$conf_file does not exist, compression not performed"; return; }

	local variables=( $(list_variables $filename) )
	local compress_args=""

	# read config file line by line
	local varname sign ndigits
	while read varname sign ndigits ; do
		# checking if name is contained in variable array with delimiting spaces can be considered safe since
		# variable do not contain spaces and if for some reason they did, the script would have broke long before :)
		if [[ $sign == "=>" ]] && [[ " ${variables[@]} " =~ " $varname " ]] ; then
			compress_args="${compress_args} --ppc ${varname}=.$ndigits"
		fi
	done <<< $(grep -v "^\ *#\|^\ *$" "$conf_file")

	# if there is something to compress, compress it
	if [[ ! -z ${compress_args} ]] ; then
		log "$LOG_FILE" "$filename | args:$compress_args"
			$NCKS -O $compress_args "${filename}" "${output_file}" &>> $LOG_FILE ||
			{ log_error "$LOG_FILE" "Failed to compress variables in $filename" && return 1; }
	else
		# copy the file for the sake of completion since this function can be used to move the output file into the right place
		[[ ${filename} != ${output_file} ]] && cp ${filename} ${output_file}
		log "$LOG_FILE" "No variables in $conf_file match. Compress not performed."
		return
	fi
}
