#!/bin/bash

# wfetchbuoy_copernicus.sh
#
# VERSION AND LAST UPDATE:
# v1.0  04/04/2022
# v1.1  01/20/2023
# v1.2  01/25/2023
# v1.3  10/01/2024
#
# PURPOSE:
#  Download Copernicus buoy database. See wget lines
#  Copernicus data
#   https://marine.copernicus.eu/
#
# USAGE:
#  Examples (from linux/terminal command line):
#   nohup sh wfetchbuoy_copernicus.sh >> nohup_wfetchbuoy_copernicus.out 2>&1 &
#
# OUTPUT:
#  Time-Series and Spectral data (netcdf format)
#
# DEPENDENCIES:
#  wget
#
# AUTHOR and DATE:
#  04/04/2022: Ricardo M. Campos, first version.
#  01/20/2023: Ricardo M. Campos, ftp address has changed. See:
# https://marine.copernicus.eu/media/pdf/November2022_Transition_Document.pdf/open
# https://data.marine.copernicus.eu/product/INSITU_GLO_WAV_DISCRETE_MY_013_045/description
# https://doi.org/10.17882/70345
# https://catalogue.marine.copernicus.eu/documents/PUM/CMEMS-INS-PUM-013-045.pdf
# https://help.marine.copernicus.eu/en/articles/4521873-how-to-download-a-dataset-from-ftp-server#h_5a34fdd0c5
#
#  01/25/2023: Ricardo M. Campos, get_buoydata_copernicus.sh renamed to wfetchbuoy_copernicus.sh
#  10/01/2024: Ricardo M. Campos, replacing the old ftp address as it does not work anymore.
# 
# PERSON OF CONTACT:
#  Ricardo M Campos: ricardo.campos@noaa.gov
#

# Time: series GL_TS_*
# Spectrum GL_WS_*

# You can download it manually via browser:
# https://data.marine.copernicus.eu/product/INSITU_GLO_WAV_DISCRETE_MY_013_045/files?subdataset=cmems_obs-ins_glo_wav_my_na_PT1H_202311--ext--history

wget "https://www.seanoe.org/data/00591/70345/data/112711.tar.gz"
wait $!
tar -zxvf 112711.tar.gz

