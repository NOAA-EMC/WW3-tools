#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
extract_Scatterometer.py

VERSION AND LAST UPDATE:
 v1.0  05/16/2022

PURPOSE:
 Script to take scatterometer swaths and collocate into a moving
  track lat/lon/time array (for ex. drifting buoy, ship etc) 
 A total of 7 satellite missions are listed below. The period of each
  scatterometer can be verified at:
  https://doi.org/10.1175/JTECH-D-19-0119.1

USAGE:
 This program processes Scatterometer data from AODN and collocate into given
  lon/lat/time array, enter as input txt file.
 Scatterometers must have been previously downloaded (see wfetchsatellite_AODN_Scatterometer.sh)
 Path where scatterometer data is saved must be informed and 
  edited (see dirs below)
 Check the pre-selected parameters below for the scatterometer collocation 
  and date interval (datemin and datemax)
 Example (from linux terminal command line):
   nohup python3 extract_scatterometer.py >> nohup_extract_scatterometer.out 2>&1 &

OUTPUT:
 Text file output_scatterometer.txt containing the collocated scatterometer
  data into lon/lat/time grid points given.
 U10: 10-meter wind speed.

DEPENDENCIES:
 See setup.py and the imports below.
 AODN scatterometer data previously downloaded (see wfetchsatellite_AODN_Scatterometer.sh)

AUTHOR and DATE:
 05/09/2022: Ricardo M. Campos, first version.

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import math
import numpy as np
import statistics
import netCDF4 as nc
import pandas as pd
import time
from calendar import timegm
import datetime as dt
import sys
import warnings; warnings.filterwarnings("ignore")
# netcdf format
fnetcdf="NETCDF4"

# Input text file name (with path if not in the same local dir)
fname='interp_AllBremChic_to_XBT.dat'  # containing lon/lat/time to be extracted

# Directory where AODN scatterometer data is saved, downloaded using wfetchsatellite_AODN_Scatterometer.sh
# dirs='/work/noaa/marine/ricardo.campos/data/AODN/scatterometer'
dirs='/media/data/observations/satellite/scatterometer/AODN_scat'
# Minimum distance (km) from the coast
mindfc=30. # in Km
# Maximum distance (m) for the collocation average
dlim=50000. # in m
# Maximum temporal distance (s) for the collocation average
maxti=2700.
# power of initial array 10**pia (size) that will be used to allocate satellite data (faster than append)
pia=10
# Satellite missions available at AODN dataset.
sdname=np.array(['ERS1','ERS2','METOPA','METOPB','OCEANSAT2','QUIKSCAT','RAPIDSCAT'])
sname=np.array(['ERS-1','ERS-2','METOP-A','METOP-B','OCEANSAT-2','QUIKSCAT','RAPIDSCAT'])
# Scatterometer Quality Control parameters
wspmax=80.


# --- INPUT ARRAY Information ---
# read with Pandas,lon/lat/time
ds = pd.read_csv(fname,comment='#',delimiter=r"\s+")
glat = np.array(ds.values[:,1]).astype('float')
glon = np.array(ds.values[:,0]).astype('float')
# Column of date strings
df = pd.DataFrame({'date_str': ds.values[:,2]})
# Define a function to convert a date string to a Unix timestamp in UTC
def date_str_to_unix_utc(date_str):
    dt_obj = dt.datetime.strptime(date_str, '%Y%m%dT%H%M%S')
    unix_timestamp_utc = int(dt_obj.replace(tzinfo=dt.timezone.utc).timestamp())
    return unix_timestamp_utc

# Apply the function to the entire column
gtime = np.array(df['date_str'].apply(date_str_to_unix_utc)).astype('double'); del df
# glon, glat, and gtime hold the information that will be used to extract the satellite data (when/where it is available)


# --- SATELLITE DATA ---
# Distance (km) between two points. A function that will be used during collocation loop.
def distance(lat1, lon1, lat2, lon2):
    # convert decimal degrees to radians
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    # haversine formula
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
    c = 2 * math.asin(math.sqrt(a))
    r = 6371  # radius of Earth in kilometers
    return c * r

# Sat files (squares) considering the lat lon of interest, for the AODN file names
auxlat=np.array(np.arange(np.floor(np.nanmin(glat))-1.,np.ceil(np.nanmax(glat))+1.,1)).astype('int')
auxlon=np.array(np.arange(0.,360.+1,1)).astype('int'); auxlon[auxlon<0]=auxlon[auxlon<0]+360

# Read and allocate satellite data into arrays
ast=np.double(np.zeros((10**pia),'d')); adfc=np.zeros((10**pia),'f'); asatid=np.zeros((10**pia),'f')*np.nan
aslat=np.zeros((10**pia),'f'); aslon=np.zeros((10**pia),'f')
awndcal=np.zeros((10**pia),'f'); awnddir=np.zeros((10**pia),'f')

ii=0
for s in range(0,np.size(sdname)):
	for j in auxlat:
		for k in auxlon:

			if j>=0:
				hem='N'
			else:
				hem='S'

			try: 
				fu=nc.Dataset(dirs+'/'+sdname[s]+'/IMOS_SRS-Surface-Waves_M_Wind-'+sname[s]+'_FV02_'+str(np.abs(j)).zfill(3)+hem+'-'+str(k).zfill(3)+'E-DM00.nc')
			except:
				print(dirs+'/'+sdname[s]+'/IMOS_SRS-Surface-Waves_M_Wind-'+sname[s]+'_FV02_'+str(np.abs(j)).zfill(3)+hem+'-'+str(k).zfill(3)+'E-DM00.nc does not exist'); vai=0
			else:
				st=np.double(fu.variables['TIME'][:]*24.*3600.+float(timegm( time.strptime('1985010100', '%Y%m%d%H') )))
				if (np.size(st)>10) and (np.nanmax(st)>np.nanmin(gtime)) and (np.nanmin(st)<np.nanmax(gtime)):
					slat=fu.variables['LATITUDE'][:,:]; sarr=int(slat.shape[0]*slat.shape[1])
					slon=fu.variables['LONGITUDE'][:,:]
					wndcal=fu.variables['WSPD_CAL'][:,:]
					swndqc=fu.variables['WSPD_CAL_quality_control'][:,:]
					try:
						wnddir=fu.variables['WND_DIR'][:]					
					except:
						print(' error reading WND_DIR')
						wnddir=wndcal*0.+999.

					# reshape and Exclude poor quality controled data
					swndqc=np.reshape(swndqc,(sarr))
					st2=np.double(np.copy(np.atleast_2d(st).T))
					for rr in range(0,slat.shape[1]-1):
						st2=np.double(np.c_[st2,np.atleast_2d(st).T])

					del st; st=np.double(np.copy(st2)); del st2
					st=np.array(np.reshape(st,(sarr))[swndqc==1]).astype('double')
					slat=np.array(np.reshape(slat,(sarr))[swndqc==1]); slon=np.array(np.reshape(slon,(sarr))[swndqc==1])
					wndcal=np.array(np.reshape(wndcal,(sarr))[swndqc==1])
					wnddir=np.array(np.reshape(wnddir,(sarr))[swndqc==1])
					del swndqc
					#
					if ii+np.size(st) <= ast.shape[0] :
						if ( (st.shape[0]==wndcal.shape[0]) & (slat.shape[0]==slon.shape[0]) ) :
							ast[ii:ii+st.shape[0]]=np.array(st).astype('double')
							aslat[ii:ii+st.shape[0]]=np.array(slat).astype('float'); aslon[ii:ii+st.shape[0]]=np.array(slon).astype('float')
							awndcal[ii:ii+st.shape[0]]=np.array(wndcal).astype('float')
							awnddir[ii:ii+st.shape[0]]=np.array(wnddir).astype('float'); 
							asatid[ii:ii+st.shape[0]]=np.array(np.zeros(wndcal.shape[0],'f')+ np.array(s)).astype('float')
							ii=ii+st.shape[0]
					else:
						sys.exit('Small array to allocate the satellite data! Increase the power of initial array (pia)')

					del st,slat,slon,wndcal,wnddir
					fu.close(); del fu

					# print(' auxLat '+repr(j)+' /  auxLon '+repr(k))
					# print(' Data array np.size '+repr(ii))
					# print(' ')

	print(' Done reading and allocating satellite data '+sdname[s])

del ii

# Simplified Quality Control Check ----
adatemin=np.nanmin(gtime)-3600.; adatemax=np.nanmax(gtime)+3600.
indq = np.where( (awndcal>0.2) & (awndcal<wspmax) & (ast>=adatemin) & (ast<=adatemax) )     

if np.size(indq)>10:
	ii=0 # apply quality control
	ast=np.double(np.copy(ast[indq[0]])); asatid=np.copy(asatid[indq[0]])
	aslat=np.copy(aslat[indq[0]]); aslon=np.copy(aslon[indq[0]])
	aslon[aslon>180.]=aslon[aslon>180.]-360.
	awndcal=np.copy(awndcal[indq[0]]); awnddir=np.copy(awnddir[indq[0]])
	# initialize arrays
	fsatid=np.zeros((gtime.shape[0]),'f')*np.nan
	fwndcal=np.zeros((gtime.shape[0]),'f')*np.nan; fwnddir=np.zeros((gtime.shape[0]),'f')*np.nan

	for t in range(0,gtime.shape[0]):
		indt = np.where( abs(ast[:]-gtime[t]) < maxti )
		if np.size(indt)>2:

			# distance
			sdist=np.zeros((np.size(indt)),'f')*np.nan
			for i in range(0,np.size(indt)):
				sdist[i]=distance(aslat[indt][i],aslon[indt][i],glat[t],glon[t])

			# within search radius
			inds=np.where(sdist<=dlim/1000.)
			if np.size(inds)>=2:
				fwndcal[t]=np.float(np.nanmean(awndcal[indt][inds]))
				fwnddir[t]=np.float(np.nanmean(awnddir[indt][inds]))
				fsatid[t]=int(statistics.mode(asatid[indt][inds]))

			del inds,sdist

		del indt
		print('  -- collocation '+repr(t)+' . '+repr(gtime.shape[0]))

	del ast, aslat, aslon, awndcal, awnddir, asatid, indq

	# Final quality control (double check)
	ind=np.where( (fwndcal<0.01) | (fwndcal>wspmax) )
	if np.any(ind):
		fwndcal[ind[0]]=np.nan
		fwnddir[ind[0]]=np.nan 
		del ind

	ind=np.where( (fwnddir<-180.) | (fwnddir>360.) )
	if np.any(ind):
		fwnddir[ind[0]]=np.nan; del ind

	print(' Total scatterometer collocated records: '+str(np.size(np.where(fwndcal>0.01)))); print(' ')
	# add header
	ds.columns = ['Lon', 'Lat', 'Time'] + list(ds.columns[3:])
	# include results in the pandas datafram
	fsatid[np.isnan(fsatid)==True]=999.; fsatid=np.array(fsatid).astype('int')
	fsatid=np.array(fsatid).astype('str'); fsatid_org=np.copy(fsatid)
	for s in range(0,np.size(np.unique(fsatid_org))-1):
		ind=np.where(fsatid_org==np.unique(fsatid_org)[s])
		if np.size(ind)>0:
			fsatid[ind[0]]=sdname[int(np.unique(fsatid_org)[s])]

	fsatid[fsatid=='999']='NaN'; del fsatid_org		
	ds['Sat'] = fsatid
	ds['U10'] = fwndcal; ds['Dir'] = fwnddir
	# save output file
	ds.to_csv('output_scatterometer.txt', sep='\t', header=True, float_format='%.3f', index=False, na_rep='NaN')

	print(' ')
	print('Output file ok and saved.')

	del ds

else:
	sys.exit(' No satellite records within the given time/date range and quality control parameters selected.')


