#!/bin/bash

# get_ww3spec_c00.sh
#
# VERSION AND LAST UPDATE:
# v1.0  04/04/2022
#
# PURPOSE:
# Operational Download and Plot of NCEP/NOAA wave spectra using python
#
# USAGE:
#  Mandatory Inputs: work directory and Spectral pointID 
#  The NCEP station/point IDs can be found at:
#   https://polar.ncep.noaa.gov/waves/viewer.shtml?-gfswave-
#  Users must have an active python. If it is not, it can be activated 
#   hrough the the souce command below (uncoment and edit line 35)
#  The python code ww3pointspec.py is used for the plots.
#  Examples (from linux/terminal command line):
#   ./get_ww3spec_c00.sh /home/user/work 41002
#   nohup ./get_ww3spec_c00.sh /home/user/work 41002 >> nohup_get_ww3spec_c00.out 2>&1 &
#
# OUTPUT:
#  A ww3spec_${ANO}${MES}${DIA}${HORA} directory is created inside the 
#   given work directory, where the png figures containing the spectrum
#   for the specific point are saved 
#  If you want to change the resolution (for publications), edit savefig
#
# DEPENDENCIES:
#  wget, tar, and python3
#  the python code ww3pointspec.py is called and must be in the same
#  directory (or path informed, or symbolic link)
#
# AUTHOR and DATE:
#  04/04/2022: Ricardo M. Campos, first version.
#
# PERSON OF CONTACT:
#  Ricardo M Campos: ricardo.campos@noaa.gov
#

source /etc/bash.bashrc

# start python
# source /home/user/python/anaconda/setanaconda3.sh

# server address 
SERVER="https://www.ftp.ncep.noaa.gov/data/nccf/com/gfs/prod"
# work directory
# DIR="/home/user/work"
DIR="$1"
# Spectral point of interest, see https://polar.ncep.noaa.gov/waves/viewer.shtml?-gfswave-
# spt="41002"
spt="$2"

# initial date cycle for the ftp
ANO=`date +%Y`
MES=`date +%m`
DIA=`date +%d`
HORA="00" # first cycle 00Z

cd ${DIR}
# Download NCEP/NOAA sprecta
wget --no-check-certificate --no-proxy -l1 -H -t1 -nd -N -np -erobots=off --tries=3 "${SERVER}/gfs."${ANO}${MES}${DIA}"/"${HORA}"/wave/station/gfswave.t00z.spec_tar.gz" -O "${DIR}/ww3spec_${ANO}${MES}${DIA}${HORA}.tgz" >> ${DIR}/logWW3spec_${ANO}${MES}${DIA}${HORA} 2>&1
#
tar -zxvf ${DIR}"/ww3spec_"${ANO}${MES}${DIA}${HORA}".tgz" >> ${DIR}/logWW3spec_${ANO}${MES}${DIA}${HORA} 2>&1
# Run python to plot spec
python3 ww3pointspec.py "${DIR}/gfswave."${spt}".spec" ${spt} >> ${DIR}/logWW3spec_${ANO}${MES}${DIA}${HORA} 2>&1

# Final organizing
mkdir "${DIR}/ww3spec_"${ANO}${MES}${DIA}${HORA} >> ${DIR}/logWW3spec_${ANO}${MES}${DIA}${HORA} 2>&1
mv "${DIR}/gfswave.${spt}.spec" "${DIR}/ww3spec_"${ANO}${MES}${DIA}${HORA} >> ${DIR}/logWW3spec_$ANO$MES$DIA$HORA 2>&1
mv "${DIR}/wspectrum_${spt}*.png" "${DIR}/ww3spec_"${ANO}${MES}${DIA}${HORA} >> ${DIR}/logWW3spec_${ANO}${MES}${DIA}${HORA} 2>&1

