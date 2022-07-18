#!/bin/bash

# get_AODN_ScatData.sh
#
# VERSION AND LAST UPDATE:
# v1.0  04/04/2022
#
# PURPOSE:
#  Script to download AODN scatterometer data. See the available files at  
#  http://thredds.aodn.org.au/thredds/catalog/IMOS/SRS/Surface-Waves/Wind-Scatterometry-DM00/catalog.html
#  Scatterometers:
#  OCEANSAT-2 ERS-1 ERS-2 METOP-A METOP-B QUIKSCAT RAPIDSCAT
#  Satellite data from Integrated Marine Observing System (IMOS), Australian Ocean Data Network (AODN)
#  https://portal.aodn.org.au/
#  Altimeter
#  https://doi.org/10.1038/s41597-019-0083-9
#  Scatterometer
#  https://doi.org/10.1175/JTECH-D-19-0119.1
#
# USAGE:
#  Three arguments are read:
#   (1) scatterometer name (exact names, see above)
#   (2) destination path
#   (3) Hemisphere, N or S
#  So if you want to download the whole database, you have to run
#   this code for each satellite and hemisphere.
#  Examples (from linux/terminal command line):
#   nohup ./get_AODN_ScatData.sh RAPIDSCAT /media/data/observations/satellite/scatterometer/AODN_scat/RAPIDSCAT S >> nohup_RAPIDSCAT_HS.out 2>&1 &
#
# OUTPUT:
#  multiple AODN satellite data (netcdf format) saved in the given 
#   directory.
#
# DEPENDENCIES:
#  wget
#
# AUTHOR and DATE:
#  04/04/2022: Ricardo M. Campos, first version.
#
# PERSON OF CONTACT:
#  Ricardo M Campos: ricardo.campos@noaa.gov
#

fname=http://thredds.aodn.org.au/thredds/fileServer/IMOS/SRS/Surface-Waves/Wind-Scatterometry-DM00
DIR="$2"

s="$1"; h="$3"
for lon in `seq -f "%03g" 0 20 340`; do
  for lat in `seq -f "%03g" 0 20 80`; do
    for adlat in `seq -f "%03g" 0 20`; do
      if [[ "$h" == "N" ]];then
        lat2=(`expr $lat + $adlat`)
      else
        lat2=(`expr $lat - $adlat`)
      fi
      for adlon in `seq -f "%03g" 0 19`; do
        lon2=(`expr $lon + $adlon`)
        test -f $DIR/IMOS_SRS-Surface-Waves_M_Wind-${s}_FV02_"$(printf "%03d" ${lat2/#-})"${h}-"$(printf "%03d" $lon2)"E-DM00.nc
        TE=$?
        if [ "$TE" -eq 1 ]; then
          wget -l1 -H -t1 -nd -N -np -erobots=off --tries=3 $fname/${s}/${lat}${h}_${lon}E/IMOS_SRS-Surface-Waves_M_Wind-${s}_FV02_"$(printf "%03d" ${lat2/#-})"${h}-"$(printf "%03d" $lon2)"E-DM00.nc -O $DIR/IMOS_SRS-Surface-Waves_M_Wind-${s}_FV02_"$(printf "%03d" ${lat2/#-})"${h}-"$(printf "%03d" $lon2)"E-DM00.nc
          wait $!
          sleep 1
          echo IMOS_SRS-Surface-Waves_M_Wind-${s}_FV02_"$(printf "%03d" ${lat2/#-})"${h}-"$(printf "%03d" $lon2)"E-DM00.nc >> listDownloaded_${s}.txt
          find $DIR -empty -type f -delete
        fi
      done
    done
  done
done
echo " " >> listDownloaded_${s}.txt
echo " OK, "${s}"   " >> listDownloaded_${s}.txt
echo " " >> listDownloaded_${s}.txt

find $DIR -empty -type f -delete

