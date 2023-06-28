#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
prep_ww3tools.py

VERSION AND LAST UPDATE:
 v1.1  01/26/2023

PURPOSE:
 This script check the ww3tools installation, set paths, download observations (minimum NDBC), and run regtests.
 prep_ww3tools.py consists of a final installation step of ww3tools.

USAGE:
 python3 prep_ww3tools.py

OUTPUT:
 Installation and functionalities of ww3tools will be tested,
  the directories and data for grid utilities, auxiliar files, and observations
  will be created and downloaded.

DEPENDENCIES:
 See setup.py and the imports below.

AUTHOR and DATE:
 01/26/2023: Ricardo M. Campos

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import os
import sys
import requests
import urllib
import numpy as np
import shutil

print(" ")
print(" ===================== WW3TOOLS Setup and preparation ====================")
print(" ")

# Check Python version  -----------------
if sys.version_info.major>=3:
	print(' OK: python version.')
else:
	sys.exit("Python Version >= 3 is required.")

# Check Installation --------------------
try:
	import ww3tools
except:
	sys.exit("Cannot find ww3tools. Please check installation at https://github.com/NOAA-EMC/WW3-tools")
	# pip install .
else:
	print(' OK: ww3tools is installed.')

try:
	from ww3tools import wread, mvalstats, pvalstats
except:
	sys.exit("Cannot find wread,mvalstats,pvalstats. Please check installation at https://github.com/NOAA-EMC/WW3-tools")
	# pip install .
else:
	print(' OK: ww3tools main functions are successfully installed.')			

# --------------------------------------

# Check internet -----------------------
try:
	requests.head('http://www.google.com/', timeout=5)
	urllib.request.urlopen('http://google.com')
except requests.ConnectionError:
	sys.exit("No internet connection available.")
else:
	print(' OK: internet connection available.')			

# --------------------

# If user do not want to use ww3_tools for ww3 validation or any measurements plot and analysis, can skip this part.
print(" ")
print(" WW3-tools can be used for simple WW3 visualization, or additional tools can be included \
for observation plots (buoys and satellie) and WW3 validation. \
For this, users need to setup the observation paths and download the measurements.")

odatareply = str(input("Do you want to continue to setup the observation data? [yes,no] :  ")).capitalize()
if 'N' in odatareply:
	print(" ")
	sys.exit("WW3-tools prep concluded.")

# Organize directories and paths for Observations --------------------
print(" ")
ghpath=str(ww3tools.__path__[0])+'/' # github path
print(" WW3-tools repository path (where you git cloned it): "+ghpath)
pspath=str(ww3tools.__path__[-1])+'/' # python system path after installation
print(" WW3-tools python installation path: "+pspath)

print(' ')
opath = str(input("Please enter the full path where new data (observations) will be saved:  "))
if len(opath)>0:
	if opath[-1]!='/':
		opath=opath+'/'
else:
	sys.exit("Observations Path must be entered.")

print(" "); print(" Path "+opath+" is where the observations will be downloaded and organized.")

wtpaths=np.array(['gridinfo','gridinfo/shapefiles','gridinfo/shapefiles/GlobalOceansSeas',
	'gridinfo/shapefiles/NOAA/HighSeasMarineZones','gridinfo/shapefiles/NOAA/OffshoreMarineZones',
	'buoys','satellite','buoys/NDBC','buoys/Copernicus','buoys/NDBC/wparam',
	'buoys/NDBC/spec','buoys/Copernicus/wparam','buoys/Copernicus/spec','satellite/AODN'])

try:
	for i in range(0,wtpaths.shape[0]):
		if not os.path.exists(opath+wtpaths[i]):
			os.makedirs(opath+wtpaths[i])
		else:
			print('  Warning: Path '+opath+wtpaths[i]+' already created. Please check the content as new data will be saved.')
except:
	sys.exit("Directories could not be created. Please check permissions.")
else:
	print(" ");print(' OK: Directories for observations created and organized.')

from ww3tools.downloadobs import ww3toolsobs, wfetchbuoy
# Include path
f = open(pspath+'downloadobs/ww3toolsobs.py','a')
f.write('\n')
f.write("def path():"); f.write('\n')
f.write("	global obsww3tools"); f.write('\n')
f.write("	obsww3tools='"+opath+"'"); f.write('\n')
f.write("	return obsww3tools"); f.write('\n')
f.writelines('\n')
f.close(); del f
print(' OK: Paths for observations included. This can be updated in the future by editing the file: '+pspath+'downloadobs/ww3toolsobs.py')
# from ww3tools import ww3toolsobs
# opath=ww3toolsobs.path()
# --------------------

# DOWNLOAD METOCEAN DATA --------------------
print(' ')
print(' ========================= DOWNLOAD OBSERVATIONS =========================')
print(' 1) GrindInfo: Bathymetry, Distance to Coast, Shapefiles of Ocean Names and Forecast Areas, and cyclone tracks (IBTrACS).')

# GrindInfo: Bathymetry, Distance to Coast, Forecast Area Shapefiles
from urllib import request
gifiles=['gridinfo/distFromCoast.nc','gridinfo/etopo1.nc','gridinfo/etopo1.nc',
	'gridinfo/shapefiles/GlobalOceansSeas/goas_v01.cpg',
	'gridinfo/shapefiles/GlobalOceansSeas/goas_v01.dbf',
	'gridinfo/shapefiles/GlobalOceansSeas/goas_v01.prj',
	'gridinfo/shapefiles/GlobalOceansSeas/goas_v01.shp',
	'gridinfo/shapefiles/GlobalOceansSeas/goas_v01.shx',
	'gridinfo/shapefiles/NOAA/HighSeasMarineZones/hz30jn17.dbf',
	'gridinfo/shapefiles/NOAA/HighSeasMarineZones/hz30jn17.prj',
	'gridinfo/shapefiles/NOAA/HighSeasMarineZones/hz30jn17.shp',
	'gridinfo/shapefiles/NOAA/HighSeasMarineZones/hz30jn17.shx',
	'gridinfo/shapefiles/NOAA/OffshoreMarineZones/oz22mr22.CPG',
	'gridinfo/shapefiles/NOAA/OffshoreMarineZones/oz22mr22.dbf',
	'gridinfo/shapefiles/NOAA/OffshoreMarineZones/oz22mr22.prj',
	'gridinfo/shapefiles/NOAA/OffshoreMarineZones/oz22mr22.shp',
	'gridinfo/shapefiles/NOAA/OffshoreMarineZones/oz22mr22.shx']

for i in range(0,len(gifiles)):
	request.urlretrieve('http://noaa-nws-gefswaves-reforecast-pds.s3.amazonaws.com/GEFSv12/'+gifiles[i],opath+gifiles[i])

# Cyclone Tracks
request.urlretrieve('https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r00/access/csv/ibtracs.since1980.list.v04r00.csv',opath+wtpaths[0]+'/ibtracs.since1980.list.v04r00.csv')

print(' 1) OK GrindInfo: successfully completed.'); print(' ')
# -----

# Buoys
print(' 2) NDBC metocean buoys.')
shutil.copyfile(ghpath+'downloadobs/allbstations.txt',opath+'buoys/allbstations.txt')
input("A list of buoys to download is given at "+opath+"buoys/allbstations.txt. \
 Feel free to keep or edit this list by removing/including stations. [ok] ")
iy = input("Initial year for download:  ")
fy = input("Final year for download:  ")

#    -- NDBC -- 
print(' '); print(' Downloading NDBC buoy data (this may take a while) ...')
wfetchbuoy.ndbc_nc(int(iy),int(fy),opath+'buoys/allbstations.txt',opath+'buoys/NDBC/wparam',opath+'buoys/NDBC/spec')

print(' 2) OK NDBC buoy data successfully downloaded at '+opath+'buoys/NDBC/')

#    -- Copernicus -- 
print(" ")
print(" The next buoy dataset, Copernicus, requires a registration and user/password , see: \
 https://data.marine.copernicus.eu/register \
 https://data.marine.copernicus.eu/product/INSITU_GLO_WAV_DISCRETE_MY_013_045/description")

creply = str(input("Do you have a data.marine.copernicus user and password ? [yes,no] :  ")).capitalize()
if 'N' in creply:
	print(" ")
	print(" Skipping Copernicus buoys ..."); print(" ")
else:
	cuser = str(input("data.marine.copernicus username:  "))
	cpassword = str(input("data.marine.copernicus password:  "))
	input("Again, feel free to edit the list of buoys at "+opath+"buoys/allbstations.txt \
	 or you can even download the entire Coperdinus buoy dataset with ww3tools/downloadobs/wfetchbuoy_copernicus.sh. [ok]")
	print(' Downloading Copernicus buoy data (this may take a while) ...')
	wfetchbuoy.copernicus_tseriesnc(cuser,cpassword,opath+'buoys/allbstations.txt',opath+'buoys/Copernicus/wparam')
	wfetchbuoy.copernicus_specnc(cuser,cpassword,opath+'buoys/allbstations.txt',opath+'buoys/Copernicus/spec')
	print(' 3) OK Copernicus buoy data successfully downloaded at '+opath+'buoys/Copernicus/'); print(' ')

# Altimeter Data
print(' 4) Satellite data requires a long time to download and therefore is not processed in this prep script. \
Using AODN quality controlled and calibrated dataset is recommended. \
Scripts to download it are wfetchsatellite_AODN_Altimeter.sh and wfetchsatellite_AODN_Scatterometer.sh, which can be run \
multiple times for each satellite mission (and hemisphere, N and S) to speed up the download. Total download time depends \
on the internet speed and how many missions it will be downloaded, varying from hours to a few days. Storage space can go \
up to 450G for Scatterometers and 120G for Altimeters. Please save files at '+opath+wtpaths[-1]); print(' ')

# Run regtests to confirm everything is working fine

print(" ")
print(" ===================== End of WW3TOOLS Setup and preparation ====================")
print(" ")

