#!/usr/bin/env python
# coding: utf-8

#### -*- coding: utf-8 -*-
'''
# Created By  : Jake Campbell
# Created Date: Thu June 17 2021

-------------------------------------------------------------
 This program downloads the historical NDBC data for a user
 provided time period. 
 -------------------------------------------------------------
 
 INPUT:      

 Directories: Observation files save locations and filenames('.nc') 
 and where to save NetCDF wave or spectra datasets. 
 
 Station name text list file locations
 
 Year, start, and end times to pull data from a specific year
 as well as specific time period inside that year in order to
 save storage space
 
 Year, Start time, End time, Interval: string date time format:
 i.e. Year: 2021 for data in the year 2021
      Start and End times: '2021-01-01' for January 1st, 2021
      Interval: '10T' for every 10 min data, '1H' for every hour data

 OUTPUT:
     
 Formatted and saved netcdf files for all available buoy stations
 for both wave and spectra data.
  
 -------------------------------------------------------------

#      Jake Campbell August 2021 jacob.campbell@noaa.gov     #

#      Ali Abdolali  August 2021 ali.abdolali@noaa.gov       #

-------------------------------------------------------------
'''
import os
import sys
import netCDF4
from pylab import *
import io
import urllib
import xarray as xr
import time
from tqdm import *
import numpy as np

### USER INPUT ################################################################
###############################################################################
# Save locations for both spectra and wave files
saveSpectra = '/work/noaa/marine/ricardo.campos/data/buoys/NDBC/ncformat/'
saveWave = '/work/noaa/marine/ricardo.campos/data/buoys/NDBC/ncformat/'

# Basing station names off current list of NDBC bouys, reading them from text file
f = open('/work/noaa/marine/ricardo.campos/data/buoys/NDBC/ncformat/allbstations.txt', "r")
stationList = f.read().splitlines()
      
# Change year to selected year or last 45 days of data (i.e. '9999' for last 45 days)
# year = input("Please enter a year to pull data from (format: YYYY):\n")
year = np.int(sys.argv[1])

# To reduce file size, select start and end date to filter data by
# start = input("Please enter start date (format: 'YYYY-MM-DD'):\n")
# end = input("Please enter end date (format: 'YYYY-MM-DD'):\n")
start = np.str(year)+'-01-01'
end = np.str(year+1)+'-01-01'

###############################################################################

# While loop to loop through Station List for both wave and spectra NetCDF 
# files and if they exist for given year (tqdm is simply a progress bar)
with tqdm(total=len(stationList), desc="{Downloading NetCDF Files}", position=0, leave=True) as pbar:

	i = 0
	for i in (range(len(stationList))):
		time.sleep(.01)
		pbar.update(1)
		stat = stationList[i]
		webWave = 'https://dods.ndbc.noaa.gov/thredds/fileServer/data/stdmet/'+str(stat)+'/'+str(stat)+'h'+str(year)+'.nc'
		webSpectra = 'https://dods.ndbc.noaa.gov/thredds/fileServer/data/swden/'+str(stat)+'/'+str(stat)+'w'+str(year)+'.nc'
		print(str(stat)+' '+str(year))
		# Set save location for each
		saveloc = saveSpectra + str(stat) +'w'+ str(year) +'.nc'
		saveloc2 = saveWave + str(stat) +'h'+ str(year) +'.nc'

		# Perform pull with check to see if spectra file exists. If it does, skip and look for next station.
		if os.path.isfile(saveloc):
			pass         
		else:
			try:
				reqSpectra = urllib.request.Request(webSpectra)
				respS = urllib.request.urlopen(reqSpectra)
			except:
				pass
			else:			
				try:
					ds_s = xr.open_dataset(io.BytesIO(respS.read()))
				except :
					os.system('wget --no-check-certificate --no-proxy -l1 -H -t1 -nd -N -np -erobots=off --tries=3 '+webSpectra)
					os.system('mv -f '+str(stat)+'w'+str(year)+'.nc '+saveSpectra)
					print('Warining with '+webSpectra+'  downloaded with wget ...')
				else:
					# Filter for specific date range before writing to NetCDF to save disk space.
					ds_s = ds_s.sel(time=slice(start, end)) 
					# Check to see if dataset is empty, if so, do not download it.
					ds_s_np = np.array(ds_s['time'])
					if ds_s_np.shape[0] == 0:
						pass
					else:
						ds_s.to_netcdf(saveloc)

                
		# Perform pull with check to see if wave file exists. If it does, skip and look for next station.    
		if os.path.isfile(saveloc2):
			pass
		else:
			try:
				reqWave = urllib.request.Request(webWave)
				respW = urllib.request.urlopen(reqWave)
			except:
				pass
			else:
				try:
					ds_w = xr.open_dataset(io.BytesIO(respW.read()))
				except :
					os.system('wget --no-check-certificate --no-proxy -l1 -H -t1 -nd -N -np -erobots=off --tries=3 '+webWave)
					os.system('mv -f '+str(stat)+'h'+str(year)+'.nc '+saveWave)
					print('Warining with '+webWave+'  downloaded with wget')
				else:
					# Can filter for specific days before writing to NetCDF.
					ds_w = ds_w.sel(time=slice(start, end))
					# Check to see if dataset is empty, if so, do not download it.
					ds_w_np = np.array(ds_w['time'])
					if ds_w_np.shape[0] == 0:
						pass
					else:
						ds_w.to_netcdf(saveloc2)

	pbar.close()
        
print('Done')

