#!/bin/bash

# download_CDIP.sh
#
# VERSION AND LAST UPDATE:
# v1.0  10/21/2024
#
# PURPOSE:
#  Script to download CDIP buoy data. See the available files at  
#  https://thredds.cdip.ucsd.edu/thredds/catalog/cdip/archive/catalog.html
#  Station info and documentation:
#  https://cdip.ucsd.edu/themes/cdip?pb=1&u2=s:209:st:1&d2=p9
#  https://docs.google.com/document/d/1Uz_xIAVD2M6WeqQQ_x7ycoM3iKENO38S4Bmn6SasHtY/edit?tab=t.0
#  https://cdip.ucsd.edu/m/products/?stn=433p1
#
# USAGE:
#  Two input arguments are required:
#   (1) A list (txt format) with the name of stations (one station per line)
#   (2) Destination path
#
#  Examples (from linux/terminal command line):
#   nohup bash download_CDIP.sh list_buoys_CDIP.txt /media/ricardo/ssdrmc/analysis/data/CDIP >> nohup_download_CDIP.out 2>&1 &
#
# OUTPUT:
#  Buoy data (netcdf format) saved in the given directory. One file per buoy ID.
#
# DEPENDENCIES:
#  wget
#
# AUTHOR and DATE:
#  10/21/2024: Ricardo M. Campos, first version.
#
# PERSON OF CONTACT:
#  Ricardo M Campos: ricardo.campos@noaa.gov
#

fname=https://thredds.cdip.ucsd.edu/thredds/fileServer/cdip/archive
# list containing buoy IDs
LIST="$1"
# output path where the netcdf files will be saved.
DIR="$2"

echo " " > $DIR/list_CDIP_Downloaded.txt
# Loop through each line (station ID) in the file
while IFS= read -r station
do

  # Download
  echo "Processing station: $station"
  wget -l1 -H -nd -N -np -erobots=off --tries=3 ${fname}/${station}p1/${station}p1_historic.nc -O $DIR/CDIP_buoy_${station}_historic.nc
  wait $!
  sleep 1
  find $DIR -empty -type f -delete
  du -sb "$DIR/CDIP_buoy_${station}_historic.nc" >> $DIR/list_CDIP_Downloaded.txt

done < "$LIST"

find $DIR -empty -type f -delete

