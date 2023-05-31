#!/bin/bash

# plot_gfsww3_forecast.sh
#
# VERSION AND LAST UPDATE:
# v1.0  05/31/2023
#
# PURPOSE:
# Download and Plot of NCEP/NOAA operational wave forecast using python.
# Wave fields of significant wave height (swh)
#
# USAGE:
#  Mandatory Inputs:
#   1) forecast cycle YYYYMMDDHH, ex. 2023053100
#   2) forecast lead time, ex. 00 (nowcast) or 24 (1-day forecast)
#   3) ww3-tools directory
#   4) work and output directory, where the grib2 forecast file and 
#     the final plot will be saved.
#
#  If you enter a forecast cycle date prior to 10-days, it will not work
#   because the ftp stores only the past 10 days of GFSWW3 operational
#   runs. Cycle hours are restricted to 00, 06, 12, and 18. The forecast
#   time is hourly until 120 hours (5 days) and then 3-hourly after that
#   until 384 hours (16 days, the limit)
#
#  Users must have an active python (or include activation in the code).
#  The python code ww3fields.py is used for the plot.
#  Example (from linux/terminal command line):
#   bash plot_gfsww3_forecast.sh 2023053106 18 /home/ricardo/WW3-tools /home/ricardo/WW3-tools/examples
#
# OUTPUT:
#  A wfields_swh_*.png plot with the global wave field of significant 
#   wave height saved at the given path. The original forecast file
#   (grib2 format) is also saved.
#
# DEPENDENCIES:
#  wget and python3
#
# AUTHOR and DATE:
#  05/31/2023: Ricardo M. Campos, first version.
#
# PERSON OF CONTACT:
#  Ricardo M Campos: ricardo.campos@noaa.gov
#

# NCEP/NOAA server address 
SERVER="https://www.ftp.ncep.noaa.gov/data/nccf/com/gfs/prod"
# forecast cycle YYYYMMDDHH, ex. 2023053000
FCYCLE="$1"
# forecast lead time
FLTIME="$2"
# ww3-tools directory, ex. /home/user/github/ww3tools
WW3TOOLSDIR="$3"
# work and output directory
WODIR="$4"
# ------------
FCYCLEDATE="${FCYCLE:0:8}"
FCYCLEHOUR="${FCYCLE:8:10}"
# output file name
OFNAME="gfswave.${FCYCLEDATE}.t${FCYCLEHOUR}z.global.0p25.f"$(printf "%03.f" ${FLTIME})".grib2"
cd ${WODIR}
# Download NCEP/NOAA forecast
wget --no-check-certificate --no-proxy -l1 -H -t1 -nd -N -np -erobots=off --tries=3 "${SERVER}/gfs.${FCYCLEDATE}/${FCYCLEHOUR}/wave/gridded/gfswave.t${FCYCLEHOUR}z.global.0p25.f"$(printf "%03.f" ${FLTIME})".grib2" -O "${WODIR}/${OFNAME}" 2>&1
# symbolic link
ln -s ${WW3TOOLSDIR}/ww3tools/ww3fields.py ${WODIR}/ww3fields.py
# Run python plot
python3 ${WODIR}/ww3fields.py ${WODIR}/${OFNAME} "swh" 2>&1
echo " "
echo "plot_gfsww3_forecast.sh completed. See "${WODIR}"/wfields_swh_*.png with the "${FLTIME}"h forecast plot for the forecast run of "${FCYCLE}
echo " "

