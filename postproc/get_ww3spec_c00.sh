#!/bin/bash

# Script to NCEP wave spectra output, and plot with python

source /etc/bash.bashrc

# start python
. /home/ricardo/python/anaconda/setanaconda3.sh

# server address
SERVER="https://www.ftp.ncep.noaa.gov/data/nccf/com/wave/prod"
# work directory
DIR=/home/ricardo/teste
# directory where figs will be saved
DIRF=/home/ricardo/teste/figs
# directory where Pleds table will be saved
DIRT=/home/ricardo/teste/pleds
# Spectral point of interest, see https://polar.ncep.noaa.gov/waves/viewer.shtml?-gfswave-
spt="SANTOS"

# initial date cycle for the ftp
ANO=`date +%Y`
MES=`date +%m`
DIA=`date +%d`
HORA="00" # first cycle 00Z

# Download NCEP/NOAA sprecta
wget --no-check-certificate --no-proxy -l1 -H -t1 -nd -N -np -erobots=off --tries=3 "${SERVER}/multi_1."${ANO}${MES}${DIA}"/multi_1.t${HORA}z.spec_tar.gz" -O $DIR/ww3spec_${ANO}${MES}${DIA}${HORA}.tgz >> $DIR/logWW3spec_$ANO$MES$DIA$HORA 2>&1
# 
tar -zxvf $DIR/ww3spec_${ANO}${MES}${DIA}${HORA}.tgz >> $DIR/logWW3spec_$ANO$MES$DIA$HORA 2>&1

# Run python to plot spec
python3 spcpointww3_textformat.py ${DIR}"/multi_1."${spt}".spec" 64 ${DIRF} ${DIRT} >> $DIR/logWW3spec_$ANO$MES$DIA$HORA 2>&1

# Final organizing
mkdir $DIR/"ww3spec_"${ANO}${MES}${DIA}${HORA} >> $DIR/logWW3spec_$ANO$MES$DIA$HORA 2>&1
mv $DIR/multi_1.*.spec $DIR/"ww3spec_"${ANO}${MES}${DIA}${HORA} >> $DIR/logWW3spec_$ANO$MES$DIA$HORA 2>&1


