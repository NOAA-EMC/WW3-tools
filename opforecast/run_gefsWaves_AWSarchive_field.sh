#!/bin/bash

# ---------------
# Retrieve AWS GEFSv12 wave forecast files from data archives
# Fetch field outputs (grib2 and converts to netcdf) for certain
#  number of days (ndays) in a MONTH and YEAR entered, using
#  gefsWaves_AWSarchive_field.sh script
# ---------------

# sh run_gefsWaves_AWSarchive_field.sh 2022 1 31 /home/ricardo/data

# module load wgrib2
# module load nco
# module load perl

YEAR="$1"
MONTH="$2"
ndays="$3" # number of days
DIRS="$4"

for DAY in $(seq $ndays); do
 cd ${DIRS}
 WTIME=${YEAR}`printf %2.2d $MONTH``printf %2.2d $DAY`
 DIRW=${DIRS}/GEFSv12waves_${WTIME}
 mkdir ${DIRW}
 ln -s ${DIRS}/get_grib.pl ${DIRW}/get_grib.pl
 ln -s ${DIRS}/get_inv.pl ${DIRW}/get_inv.pl
 sh ${DIRS}/gefsWaves_AWSarchive_field.sh ${WTIME} ${DIRW}
 wait $!
done

