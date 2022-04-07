"""
modelBuoy_collocation_forecast.py

VERSION AND LAST UPDATE:
 v1.0  04/04/2022

PURPOSE:
 Collocation/pairing ww3 point output forecast results with wave buoys.
 Matchups of ww3 results and buoy data are generated, for the same
  points (lat/lon) and time.

USAGE:
 This code is designed for ww3 forecasts, results will have two time 
  dimensions (cycle time, and forecast lead time)
 Ww3 netcdf results for point output tables (tab) is utilized.
 In order to process multiple ww3 files and append results, a list of
  ww3 output tab files is read, ww3list.txt, which must be informed. 
 It uses two public buoy databases, NDBC and Copernicus,
  which (at least one) must have been previously downloaded. See
  get_buoydata_copernicus.sh and retrieve_ndbc_nc.py
 Edit mpath, ndbcp, and copernp to set paths.
 Users have to confirm the buoys' names at the "f=nc.Dataset" lines below.
 Python code can be run directly, without input arguments.

OUTPUT:
 netcdf file ww3buoy_collocation_*.nc containing the matchups of buoy 
  and ww3 data, for the stations (lat/lon) where both data sources 
  are available.

DEPENDENCIES:
 See dependencies.py and the imports below.
 ww3 table results in netcdf format (list of files ww3list.txt)
 NDBC buoy data (see retrieve_ndbc_nc.py)
 Copernicus buoy data (see get_buoydata_copernicus.sh)

AUTHOR and DATE:
 04/04/2022: Ricardo M. Campos, first version.

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import numpy as np
from matplotlib.mlab import *
from pylab import *
import xarray as xr
import netCDF4 as nc
import time
from time import strptime
from calendar import timegm
import warnings; warnings.filterwarnings("ignore")
# netcdf format
fnetcdf="NETCDF4"

# Paths
# WW3 Model
mpath="/ww3runs/c00"
# NDBC buoys
ndbcp="/data/buoys/NDBC/wparam"
# Copernicus buoys
copernp="/data/buoys/Copernicus/wtimeseries"
# import os; os.system("ls "+mpath+"/*tab.nc > ww3list.txt &")
wlist=np.loadtxt('ww3list.txt',dtype=str)

# Read Data
# Model
for t in range(0,size(wlist)):

	try:
		ds = xr.open_dataset(wlist[t]); f=nc.Dataset(wlist[t])
	except:
		print(" Cannot open "+wlist[t])
	else:
		if t==0:
			# forecast lead time
			flt = np.array(f.variables['time'][:]*24*3600 + timegm( strptime(np.str(f.variables['time'].units).split(' ')[2][0:4]+'01010000', '%Y%m%d%H%M') )).astype('double')
			mtime=np.array([flt]).astype('double')
			sflt=np.array( flt - np.nanmin(flt) ).astype('double')
			# list of station/buoy names
			auxstationname=ds['station_name'].values[:,:]; stname=[]
			for i in range(0,auxstationname.shape[0]):
				stname=np.append(stname,"".join(np.array(auxstationname[i,:]).astype('str')))

			# Initial arrays
			# stream/test, members, point, cycle/time, forecast_lead_time
			mhs=np.zeros((size(stname),size(wlist),flt.shape[0]),'f')*np.nan
			mtm=np.zeros((size(stname),size(wlist),flt.shape[0]),'f')*np.nan
			mdm=np.zeros((size(stname),size(wlist),flt.shape[0]),'f')*np.nan
			# buoy data
			bhs=np.zeros((size(stname),size(wlist),flt.shape[0]),'f')*np.nan
			btm=np.zeros((size(stname),size(wlist),flt.shape[0]),'f')*np.nan
			bdm=np.zeros((size(stname),size(wlist),flt.shape[0]),'f')*np.nan
			lat=np.zeros(size(stname),'f')*np.nan; lon=np.zeros(size(stname),'f')*np.nan

		else:
			flt = np.array(f.variables['time'][:]*24*3600 + timegm( strptime(np.str(f.variables['time'].units).split(' ')[2][0:4]+'01010000', '%Y%m%d%H%M') )).astype('double')
			mtime = np.append(mtime,np.array([flt[0]+sflt]).astype('double'),axis=0)

		ahs = np.array(ds['hs'].values[:,:]).T
		atm = np.array(ds['tr'].values[:,:]).T
		adm = np.array(ds['th1m'].values[:,:]).T
		ds.close(); del ds
		ignore,ind,ignore = np.intersect1d(mtime[t,:], flt, return_indices=True)
		mhs[:,t,ind] = ahs
		mtm[:,t,ind] = atm
		mdm[:,t,ind] = adm
		del ahs,atm,adm

	print("done "+repr(t))

# help reading NDBC buoys, divided by year
yrange=np.array(np.arange(time.gmtime(mtime.min())[0],time.gmtime(mtime.min())[0]+1,1)).astype('int')
# Buoys
for b in range(0,size(stname)):

	ahs=[]
	try:
		ahs=[];atm=[];adm=[];atime=[]
		for y in yrange:

			f=nc.Dataset(ndbcp+"/"+stname[b]+"h"+repr(y)+".nc")
			ahs = np.append(ahs,f.variables['wave_height'][:,0,0])
			atm = np.append(atm,f.variables['average_wpd'][:,0,0])
			adm = np.append(adm,f.variables['mean_wave_dir'][:,0,0])
			atime = np.append(atime,np.array(f.variables['time'][:]).astype('double'))
			lat[b] = f.variables['latitude'][:]; lon[b] = f.variables['longitude'][:]
			f.close(); del f

	except:
		try:
			f=nc.Dataset(copernp+"/GL_TS_MO_"+stname[b]+".nc")
			ahs = np.nanmean(f.variables['VHM0'][:,:],axis=1)
			atm = np.nanmean(f.variables['VTM02'][:,:],axis=1)
			adm = np.nanmean(f.variables['VMDR'][:,:],axis=1)
			atime = np.array(f.variables['TIME'][:]*24*3600 + timegm( strptime('195001010000', '%Y%m%d%H%M') )).astype('double')
			# lat[b] = f.variables['LATITUDE'][:]; lon[b] = f.variables['LONGITUDE'][:]
			lat[b] = np.nanmean(f.variables['LATITUDE'][:]); lon[b] = np.nanmean(f.variables['LONGITUDE'][:])
			f.close(); del f
		except:
			ahs=[]

	else:
		if size(ahs>0):
			c=0
			for t in range(0,size(wlist)):
				for i in range(0,mtime.shape[1]):
					indt=np.where(np.abs(atime-mtime[t,i])<1800.)
					if size(indt)>0:
						if np.any(ahs[indt[0]].mask==False):
							bhs[b,t,i] = np.nanmean(ahs[indt[0]][ahs[indt[0]].mask==False])
							c=c+1
						if np.any(atm[indt[0]].mask==False):
							btm[b,t,i] = np.nanmean(atm[indt[0]][atm[indt[0]].mask==False])
						if np.any(adm[indt[0]].mask==False):
							bdm[b,t,i] = np.nanmean(adm[indt[0]][adm[indt[0]].mask==False])

			print("counted "+repr(c)+" at "+stname[b])

	print("done "+stname[b])
	del ahs

# Simple quality-control (range)
ind=np.where((bhs>30.)|(bhs<0.0))
if size(ind)>0:
	bhs[ind]=np.nan; del ind

ind=np.where((btm>40.)|(btm<0.0))
if size(ind)>0:
	btm[ind]=np.nan; del ind

ind=np.where((bdm>360.)|(bdm<-180.))
if size(ind)>0:
	bdm[ind]=np.nan; del ind

ind=np.where((mhs>30.)|(mhs<0.0))
if size(ind)>0:
	mhs[ind]=np.nan; del ind

ind=np.where((mtm>40.)|(mtm<0.0))
if size(ind)>0:
	mtm[ind]=np.nan; del ind

ind=np.where((mdm>360.)|(mdm<-180.))
if size(ind)>0:
	mdm[ind]=np.nan; del ind

# Save netcdf output file
lon[lon>180.]=lon[lon>180.]-360.
inidate=np.str(time.gmtime(mtime.min())[0])+np.str(time.gmtime(mtime.min())[1]).zfill(2)+np.str(time.gmtime(mtime.min())[2]).zfill(2)+np.str(time.gmtime(mtime.min())[3]).zfill(2)
findate=np.str(time.gmtime(mtime.max())[0])+np.str(time.gmtime(mtime.max())[1]).zfill(2)+np.str(time.gmtime(mtime.max())[2]).zfill(2)+np.str(time.gmtime(mtime.max())[3]).zfill(2)
# Save netcdf
ncfile = nc.Dataset("ww3buoy_collocation_"+inidate+"to"+findate+".nc", "w", format=fnetcdf) 
ncfile.history="Collocation of WW3 point output forecast (table) and NDBC/Copernicus Buoys. Total of "+repr(bhs[bhs>0.].shape[0])+" observations or pairs model/observation."
# create  dimensions. 2 Dimensions
ncfile.createDimension('buoypoints', bhs.shape[0] )
ncfile.createDimension('fcycletime', bhs.shape[1] )
ncfile.createDimension('forecastleadtime', bhs.shape[2] )
# create variables.
vt = ncfile.createVariable('time',np.dtype('float64').char,('fcycletime','forecastleadtime'))
vstname = ncfile.createVariable('buoyID',dtype('a25'),('buoypoints'))
vlat = ncfile.createVariable('latitude',np.dtype('float32').char,('buoypoints'))
vlon = ncfile.createVariable('longitude',np.dtype('float32').char,('buoypoints'))
#
vmhs = ncfile.createVariable('model_hs',np.dtype('float32').char,('buoypoints','fcycletime','forecastleadtime'))
vmtm = ncfile.createVariable('model_tm',np.dtype('float32').char,('buoypoints','fcycletime','forecastleadtime'))
vmdm = ncfile.createVariable('model_dm',np.dtype('float32').char,('buoypoints','fcycletime','forecastleadtime'))
vbhs = ncfile.createVariable('buoy_hs',np.dtype('float32').char,('buoypoints','fcycletime','forecastleadtime'))
vbtm = ncfile.createVariable('buoy_tm',np.dtype('float32').char,('buoypoints','fcycletime','forecastleadtime'))
vbdm = ncfile.createVariable('buoy_dm',np.dtype('float32').char,('buoypoints','fcycletime','forecastleadtime'))
# Assign units
vlat.units = 'degrees_north' ; vlon.units = 'degrees_east'
vt.units = 'seconds since 1970-01-01T00:00:00+00:00'
vmhs.units='m'; vbhs.units='m'
vmtm.units='s'; vbtm.units='s'
vmdm.units='degrees'; vbdm.units='degrees'
# Allocate Data
vt[:,:]=mtime[:,:]; vstname[:]=stname[:]
vlat[:] = lat[:]; vlon[:] = lon[:]
vmhs[:,:,:]=mhs[:,:,:]
vmtm[:,:,:]=mtm[:,:,:]
vmdm[:,:,:]=mdm[:,:,:]
vbhs[:,:,:]=bhs[:,:,:]
vbtm[:,:,:]=btm[:,:,:]
vbdm[:,:,:]=bdm[:,:,:]
#
ncfile.close()
print(' ')
print('netcdf ok ')

