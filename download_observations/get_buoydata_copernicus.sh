#!/bin/bash

# get_buoydata_copernicus.sh
#
# VERSION AND LAST UPDATE:
# v1.0  04/04/2022
#
# PURPOSE:
#  Download Copernicus buoy database. See wget lines
#  Copernicus data
#   https://marine.copernicus.eu/
#
# USAGE:
#  You have to enter login and password below.
#  Examples (from linux/terminal command line):
#   nohup ./get_buoydata_copernicus.sh >> nohup_get_buoydata_copernicus.out 2>&1 &
#
# OUTPUT:
#  Time-Series and Spectral data (netcdf format)
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

# Time series
fname=ftp://my.cmems-du.eu/Core/INSITU_GLO_WAVE_REP_OBSERVATIONS_013_045/history/MO/GL_TS_*
wget -c -N --user='' --password='' $fname 

# Spectrum
fname=ftp://my.cmems-du.eu/Core/INSITU_GLO_WAVE_REP_OBSERVATIONS_013_045/history/MO/GL_WS_*
wget -c -N --user='' --password='' $fname 

