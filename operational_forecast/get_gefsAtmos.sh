#!/bin/bash

########################################################################
# get_gefsAtmos.sh
#
# VERSION AND LAST UPDATE:
#   v1.0  11/09/2022
#
# PURPOSE:
#  Script to download NOAA Global Ensemble Forecast System (GEFS),
#   atmospheric operational forecast. Download from ftp, not nomads.
#
# USAGE:
#  There are two path variables to setup, one containing the scripts
#    and the other where the data will be placed. See DIRS and DIR below.
#  The variable VARSGET has the GEFS variables to download. If you 
#    change it, re-evaluate the variable TAM, associated with the 
#    expected size of each grib2 file downloaded.
#  An additional compression with ncks reads an argument definded as dp,
#    which can be modified with caution (don't use very small values).
#  By default it downloads the operational forecast of the day before,
#    defined by pa=1. If you want data for the same day, uncomment the
#    lines with YEAR, MONTH, DAY and comment the following lines with pa.
#  To reduce even more the final size of the file, a min/max latitude 
#    is defined. By default this is set from 79S to 90N. The latitudes 
#    south of 79S are excluded because there is no ocean in the area.
#  This program can be directly run with ./get_gefsWaves.sh
#
# OUTPUT:
#  One cycle of NOAA's GEFS Atmospheric operational forecast, including
#   all the 30 members plus the control member. Netcdf format.
#  A directory gefsAtmos.$YEAR$MONTH$DAY$HOUR is generated containing the
#    data, plus a log file logGEFS_$YEAR$MONTH$DAY$HOUR.
#
# DEPENDENCIES:
#  wget, perl, NCO, wgrib2
#  For most linux systems, wget and perl are already included. The
#    program NCO can be downloaded via apt-get (debian/ubuntu):
#    sudo apt-get install nco
#    the program wgrib2 can be obtained at
#    https://www.ftp.cpc.ncep.noaa.gov/wd51we/wgrib2/
#
#  Finally, the two codes get_grib.pl and get_inv.pl must be included in 
#    the same directory where you save this code get_gefsAtmos.sh
#
# AUTHOR and DATE:
#  11/09/2022: Ricardo M. Campos, first simplified version 
#
# PERSON OF CONTACT:
#  Ricardo M. Campos: ricardo.campos@noaa.gov
#
########################################################################

# Cluster module loads
# module load wgrib2
# module load nco

# path where this code and get_grib.pl and get_inv.pl are saved
DIRS=/data/archives/scripts
# path where directory will be created and downloaded files will be saved
DIR=/data/archives/data

# server address
SERVER=https://ftpprd.ncep.noaa.gov/
s1="pgrb2a" 

# convertion and compression
WGRIB2=$(which wgrib2)
NCKS=$(which ncks)
NCATTED=$(which ncatted)
NCRCAT=$(which ncrcat)

# cutoff decimals to reduce file size
dp=2
# limits for domain selection
latmin=-79.
latmax=90.

# Global Ensemble Forecast System (GEFS) https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/global-ensemble-forecast-system-gefs
# https://www.ftp.ncep.noaa.gov/data/nccf/com/gens/prod/
# variable names to be downloaded.
VARSGET=":UGRD:850 mb:|:UGRD:10 m above ground:|:VGRD:850 mb:|:VGRD:10 m above ground:|:PRMSL:mean sea level:|:HGT:850 mb:|:TMP:2 m above ground:|:RH:2 m above ground:"

# initial date cycle for the ftp
#YEAR=`date +%Y`
#MONTH=`date +%m`
#DAY=`date +%d`
pa=1
YEAR=`date --date=-$pa' day' '+%Y'`
MONTH=`date --date=-$pa' day' '+%m'`
DAY=`date --date=-$pa' day' '+%d'`
HOUR="00" # first cycle 00Z

cd $DIR
# create directory
mkdir -p $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR
# all information about fetching the grib2 files will be saved in the log file 
echo "  " > $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR

# hours of forecast lead time to be dowloaded
h1="`seq -f "%03g" 0 3 240`"
h2="`seq -f "%03g" 246 6 840`"
fleads=${h1}' '${h2}
# number of ensemble members
ensblm="`seq -f "%02g" 0 1 30`"

for h in $fleads;do
  for e in $ensblm;do

    echo "  " >> $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR
    echo " ======== GEFS Forecast: $YEAR$MONTH$DAY$HOUR  $h ========" >> $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR 

    # size TAM and tries TRIES will control the process
    TAM=0
    TRIES=1
    # while file has a lower size the expected or attemps are less than 130 (almos 11 hours trying) it does:
    while [ $TAM -lt 1000000 ] && [ $TRIES -le 130 ]; do
      # sleep 5 minutes between attemps
      if [ ${TRIES} -gt 5 ]; then
        sleep 300
      fi

      if [ ${TAM} -lt 1000000 ]; then
          echo "  " >> $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR
          echo " attempt number: $TRIES" >> $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR
          echo "  " >> $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR
          # main line where get_inv.pl and get_grib.pl are used to fech the grib2 file  
          if [ ${e} == 00 ]; then  
             $DIRS/get_inv.pl $SERVER/data/nccf/com/gens/prod/gefs.$YEAR$MONTH$DAY/$HOUR/atmos/${s1}p5/gec${e}.t${HOUR}z.${s1}.0p50.f"$(printf "%03.f" $h)".idx | egrep "($VARSGET)" | $DIRS/get_grib.pl $SERVER/data/nccf/com/gens/prod/gefs.$YEAR$MONTH$DAY/$HOUR/atmos/${s1}p5/gec${e}.t${HOUR}z.${s1}.0p50.f"$(printf "%03.f" $h)" $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/gefs${e}.t${HOUR}z.pgrb2f"$(printf "%03.f" $h)".grb2 >> $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR 2>&1         
          else         
             $DIRS/get_inv.pl $SERVER/data/nccf/com/gens/prod/gefs.$YEAR$MONTH$DAY/$HOUR/atmos/${s1}p5/gep${e}.t${HOUR}z.${s1}.0p50.f"$(printf "%03.f" $h)".idx | egrep "($VARSGET)" | $DIRS/get_grib.pl $SERVER/data/nccf/com/gens/prod/gefs.$YEAR$MONTH$DAY/$HOUR/atmos/${s1}p5/gep${e}.t${HOUR}z.${s1}.0p50.f"$(printf "%03.f" $h)" $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/gefs${e}.t${HOUR}z.pgrb2f"$(printf "%03.f" $h)".grb2 >> $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR 2>&1
          fi
          # test if the downloaded file exists
          test -f $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/gefs${e}.t${HOUR}z.pgrb2f"$(printf "%03.f" $h)".grb2 >> $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR 2>&1
          TE=$?
          if [ ${TE} -eq 1 ]; then
            TAM=0
          else
            # check size of each file
            TAM=`du -sb $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/gefs${e}.t${HOUR}z.pgrb2f"$(printf "%03.f" $h)".grb2 | awk '{ print $1 }'` >> $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR 2>&1
          fi
      fi

      TRIES=`expr $TRIES + 1`
    done

  done
done
sleep 2

cd $DIR

# Cleaning 
# --- Remove directories older than 5 days
# find gefsAtmos.?????????? -ctime +4 -type d | xargs rm -rf >> $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR 2>&1 

# Post-processing: select area and reduces resolution, to save disk space. ------------------
echo " " >> $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR
echo " Post-Processing. select area and reduces resolution " >> $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR

for h in $fleads;do
  for e in $ensblm;do

     arqn=$DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/gefs${e}.t${HOUR}z.pgrb2f"$(printf "%03.f" $h)"
     test -f ${arqn}.grb2 >> $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR 2>&1
     TE=$?
     if [ ${TE} -eq 1 ]; then
       echo " File ${arqn}.grb2 does not exist. Download failed." >> $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR
     else
       $WGRIB2 ${arqn}.grb2 -netcdf ${arqn}.saux.nc >> $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR 2>&1
       wait $!
       $NCKS -4 -L 1 -d latitude,${latmin},${latmax} ${arqn}.saux.nc ${arqn}.saux2.nc >> $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR 2>&1
       wait $!
       $NCKS --ppc default=.$dp ${arqn}.saux2.nc ${arqn}.nc >> $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR 2>&1
       wait $!      
       $NCATTED -a _FillValue,,o,f,NaN ${arqn}.nc >> $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR 2>&1
       wait $!
       rm -f ${arqn}.grb2
       rm -f ${arqn}.saux*
       echo " File ${arqn} converted to netcdf and compressed with success. " >> $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR
     fi

  done
done

# Merge all netcdf files
for e in $ensblm;do
  $NCRCAT $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/gefs${e}.t${HOUR}z.pgrb2f*.nc -O $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/gefs.$YEAR$MONTH$DAY$HOUR.m${e}.${s1}.nc >> $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/logGEFS_$YEAR$MONTH$DAY$HOUR 2>&1
  wait $!
  sleep 2
  rm -f $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR/gefs${e}.t${HOUR}z.pgrb2f*.nc
done

# permissions and groups
chmod -R 775 $DIR/gefsAtmos.$YEAR$MONTH$DAY$HOUR

