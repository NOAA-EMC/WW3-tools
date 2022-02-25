#!/bin/bash
# Download Copernicus buoy database
# You have to enter login and password below.

# Time series
fname=ftp://my.cmems-du.eu/Core/INSITU_GLO_WAVE_REP_OBSERVATIONS_013_045/history/MO/GL_TS_*
wget -c -N --user='' --password='' $fname 

# Spectrum
fname=ftp://my.cmems-du.eu/Core/INSITU_GLO_WAVE_REP_OBSERVATIONS_013_045/history/MO/GL_WS_*
wget -c -N --user='' --password='' $fname 

