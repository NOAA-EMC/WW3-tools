#!/bin/bash

########################################################################
# get_gfsWaves.sh
#
# VERSION AND LAST UPDATE:
#   v1.0  11/09/2022
#
# PURPOSE:
#  Script to download NOAA Global Forecast System (GFS), Wave 
#   Forecast from WAVEWATCH III operational. Download from ftp, not nomads.
#
# USAGE:
#  There are two path variables to setup, one containing the scripts 
#    and the other where the data will be placed. See DIRS and DIR below.
#  The variable VARSGET has the GFS variables to download. If you
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
#  This program can be directly run with ./get_gfsWaves.sh
#
# OUTPUT:
#  One cycle of NOAA's GFS Waves operational forecast. Netcdf format.
#  A directory gfsWave.$YEAR$MONTH$DAY$HOUR is generated containing the
#    data, plus a log file logGFS_$YEAR$MONTH$DAY$HOUR.
#
# DEPENDENCIES:
#  wget, perl, NCO, wgrib2
#  For most linux systems, wget and perl are already included. The
#    program NCO can be downloaded via apt-get (debian/ubuntu):
#    sudo apt-get install nco
#    the program wgrib2 can be downloaded at
#    https://www.ftp.cpc.ncep.noaa.gov/wd51we/wgrib2/
#
#  Finally, the two codes get_grib.pl and get_inv.pl must be included in 
#    the same directory where you save this code get_gfsWaves.sh
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
SERVER=https://www.ftp.ncep.noaa.gov
# s1="global.0p16"
s1="global.0p25"

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
# variable naMONTH to be downloaded.
VARSGET=":UGRD:surface:|:VGRD:surface:|:HTSGW:surface:|:PERPW:surface:|:DIRPW:surface:"

# Initial date cycle for the ftp
#YEAR=`date +%Y`
#MONTH=`date +%m`
#DAY=`date +%d`
pa=2 #  days into the past. pa=1 dowloades data from yesterday's cycle
YEAR=`date --date=-$pa' day' '+%Y'`
MONTH=`date --date=-$pa' day' '+%m'`
DAY=`date --date=-$pa' day' '+%d'`
HOUR="00" # first cycle 00Z

cd $DIR

# create directory
mkdir -p $DIR/gfsWave.$YEAR$MONTH$DAY$HOUR
# all information about fetching the grib2 files will be saved in the log file 
echo "  " > $DIR/gfsWave.$YEAR$MONTH$DAY$HOUR/logGFS_$YEAR$MONTH$DAY$HOUR

# hours of forecast to be dowloaded
h1="`seq -f "%03g" 0 1 120`"
h2="`seq -f "%03g" 123 3 384`"
fleads=${h1}' '${h2}
#
for h in $fleads;do
  echo "  " >> $DIR/gfsWave.$YEAR$MONTH$DAY$HOUR/logGFS_$YEAR$MONTH$DAY$HOUR
  echo " ======== GFS Wave Forecast: $YEAR$MONTH$DAY$HOUR  $h ========" >> $DIR/gfsWave.$YEAR$MONTH$DAY$HOUR/logGFS_$YEAR$MONTH$DAY$HOUR 
  # size TAM and tries TRIES will control the process
  TAM=0
  TRIES=1
  # while file has a lower size the expected or attemps are less than 130 (almos 11 hours trying) it does:
  while [ ${TAM} -lt 1000  ] && [ ${TRIES} -le 130 ]; do
    # sleep 5 minutes between attemps
    if [ ${TRIES} -gt 5 ]; then
      sleep 300
    fi

    if [ ${TAM} -lt 1000 ]; then
        echo "  " >> $DIR/gfsWave.$YEAR$MONTH$DAY$HOUR/logGFS_$YEAR$MONTH$DAY$HOUR
        echo " attempt number: $TRIES" >> $DIR/gfsWave.$YEAR$MONTH$DAY$HOUR/logGFS_$YEAR$MONTH$DAY$HOUR
        echo "  " >> $DIR/gfsWave.$YEAR$MONTH$DAY$HOUR/logGFS_$YEAR$MONTH$DAY$HOUR
        # main line where get_inv.pl and get_grib.pl are used to fech the grib2 file
        $DIRS/get_inv.pl $SERVER/data/nccf/com/gfs/prod/gfs.$YEAR$MONTH$DAY/$HOUR/wave/gridded/gfswave.t${HOUR}z.${s1}.f"$(printf "%03.f" $h)".grib2.idx | egrep "($VARSGET)" | $DIRS/get_grib.pl $SERVER/data/nccf/com/gfs/prod/gfs.$YEAR$MONTH$DAY/$HOUR/wave/gridded/gfswave.t${HOUR}z.${s1}.f"$(printf "%03.f" $h)".grib2 $DIR/gfsWave.$YEAR$MONTH$DAY$HOUR/gfswave.t${HOUR}z.f"$(printf "%03.f" $h)".grib2 >> $DIR/gfsWave.$YEAR$MONTH$DAY$HOUR/logGFS_$YEAR$MONTH$DAY$HOUR 2>&1         
        # test if the downloaded file exists
        test -f $DIR/gfsWave.$YEAR$MONTH$DAY$HOUR/gfswave.t${HOUR}z.f"$(printf "%03.f" $h)".grib2 >> $DIR/gfsWave.$YEAR$MONTH$DAY$HOUR/logGFS_$YEAR$MONTH$DAY$HOUR 2>&1
        TE=$?
        if [ "$TE" -eq 1 ]; then
          TAM=0
        else
          # check size of each file
          TAM=`du -sb $DIR/gfsWave.$YEAR$MONTH$DAY$HOUR/gfswave.t${HOUR}z.f"$(printf "%03.f" $h)".grib2 | awk '{ print $1 }'` >> $DIR/gfsWave.$YEAR$MONTH$DAY$HOUR/logGFS_$YEAR$MONTH$DAY$HOUR 2>&1
        fi
    fi

    TRIES=`expr $TRIES + 1`
  done
done

# Cleaning 
# --- Remove directories older than 5 days
cd $DIR
# find gfswave.?????????? -ctime +4 -type d | xargs rm -rf >> $DIR/gfsWave.$YEAR$MONTH$DAY$HOUR/logGFS_$YEAR$MONTH$DAY$HOUR 2>&1 

# Post-processing: select area, compress, and reduce decimals resolution, to save disk space. ------------------
echo " " >> $DIR/gfsWave.$YEAR$MONTH$DAY$HOUR/logGFS_$YEAR$MONTH$DAY$HOUR
echo " Post-Processing. Netcdf4 compression " >> $DIR/gfsWave.$YEAR$MONTH$DAY$HOUR/logGFS_$YEAR$MONTH$DAY$HOUR

for h in $fleads;do

  arqn=$DIR/gfsWave.$YEAR$MONTH$DAY$HOUR/gfswave.t${HOUR}z.f$h
  test -f ${arqn}.grib2 >> $DIR/gfsWave.$YEAR$MONTH$DAY$HOUR/logGFS_$YEAR$MONTH$DAY$HOUR 2>&1
  TE=$?
  if [ "$TE" -eq 1 ]; then
     echo " File ${arqn}.grib2 does not exist. Failed to download " >> $DIR/gfsWave.$YEAR$MONTH$DAY$HOUR/logGFS_$YEAR$MONTH$DAY$HOUR
  else
     $WGRIB2 ${arqn}.grib2 -netcdf ${arqn}.saux.nc >> $DIR/gfsWave.$YEAR$MONTH$DAY$HOUR/logGFS_$YEAR$MONTH$DAY$HOUR 2>&1
     wait $!
     $NCKS -4 -L 1 -d latitude,${latmin},${latmax} ${arqn}.saux.nc ${arqn}.saux2.nc >> $DIR/gfsWave.$YEAR$MONTH$DAY$HOUR/logGFS_$YEAR$MONTH$DAY$HOUR 2>&1
     wait $!
     $NCKS --ppc default=.$dp ${arqn}.saux2.nc ${arqn}.nc >> $DIR/gfsWave.$YEAR$MONTH$DAY$HOUR/logGFS_$YEAR$MONTH$DAY$HOUR 2>&1
     wait $!
     $NCATTED -a _FillValue,,o,f,NaN ${arqn}.nc >> $DIR/gfsWave.$YEAR$MONTH$DAY$HOUR/logGFS_$YEAR$MONTH$DAY$HOUR 2>&1
     wait $!
     rm -f ${arqn}.grib2
     rm ${arqn}.saux*
     echo " File ${arqn} converted to netcdf and compressed with success. " >> $DIR/gfsWave.$YEAR$MONTH$DAY$HOUR/logGFS_$YEAR$MONTH$DAY$HOUR
  fi

done

# Merge all netcdf files
$NCRCAT $DIR/gfsWave.$YEAR$MONTH$DAY$HOUR/*.nc -O $DIR/gfsWave.$YEAR$MONTH$DAY$HOUR/gfsWave.$YEAR$MONTH$DAY$HOUR.${s1}.nc
sleep 10
rm -f $DIR/gfsWave.$YEAR$MONTH$DAY$HOUR/*z.f*

# permissions and groups
chmod -R 775 $DIR/gfsWave.$YEAR$MONTH$DAY$HOUR

