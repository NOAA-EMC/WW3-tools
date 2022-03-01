# Collocation/pairing GEFS ww3 point output forecast results with wave buoys.
# This code is designed for ww3 ensemble forecasts, with three time dimensions (ensemble member, cycle time, and forecast lead time)
# two ww3 netcdf point output formats are expected, tab and spec. Only spectral point output has wind speed.

import numpy as np
from matplotlib.mlab import *
from pylab import *
import xarray as xr
import netCDF4 as nc
import sys
import os
import time
from time import strptime
from calendar import timegm
import warnings; warnings.filterwarnings("ignore")
# netcdf format
fnetcdf="NETCDF4"

# Users have to check the buoys' names at the "f=nc.Dataset" lines below.

# Paths
# WW3 Model
mpath="/work/noaa/marine/ricardo.campos/work/ww3runs/results/stream01"
# NDBC buoys
ndbcp="/work/noaa/marine/ricardo.campos/data/buoys/NDBC/ncformat/wparam"
# Copernicus buoys
copernp="/work/noaa/marine/ricardo.campos/data/buoys/Copernicus/wtimeseries"

# Ensemble members
os.system("ls "+mpath+"/ > enslist.txt")
members=np.loadtxt('enslist.txt',dtype=str)
# WW3 files inside each dir
os.system("ls -1 "+mpath+"/"+members[0]+"/*_tab.nc | xargs -n 1 basename > ww3list.txt")
wlist=np.loadtxt('ww3list.txt',dtype=str)

# READ DATA
# Model
for j in range(0,size(members)):
	for t in range(0,size(wlist)):

		try:
			# Table
			fname=mpath+"/"+members[j]+"/"+wlist[t]
			ds = xr.open_dataset(fname); f=nc.Dataset(fname)
			# Spectrum (contains wind)
			ds2 = xr.open_dataset(fname.split('_')[0]+'_spec.nc'); f2=nc.Dataset(fname.split('_')[0]+'_spec.nc')
		except:
			print(" Cannot open "+mpath+"/"+members[j]+"/"+wlist[t])
		else:

			flt2 = np.array(f2.variables['time'][:]*24*3600 + timegm( strptime(np.str(f.variables['time'].units).split(' ')[2][0:4]+'01010000', '%Y%m%d%H%M') )).astype('double')			
			flt = np.array(f.variables['time'][:]*24*3600 + timegm( strptime(np.str(f.variables['time'].units).split(' ')[2][0:4]+'01010000', '%Y%m%d%H%M') )).astype('double')

			if (flt2.min() != flt.min()) or (flt2.max() != flt.max()):
				sys.exit(' Initial/final time of tab and spec ww3 files must be the same.')
			
			if j==0 and t==0:
				
				if flt2.shape[0] > flt.shape[0]:
					auxf = np.copy(flt)
					flt = np.copy(flt2)
					flt2 = np.copy(auxf); del auxf

				# forecast lead time
				mtime=np.array([flt2]).astype('double')
				sflt=np.array( flt2 - np.nanmin(flt2) ).astype('double')
				
				# list of station/buoy names
				auxstationname=ds['station_name'].values[:,:]; stname=[]
				for i in range(0,auxstationname.shape[0]):
					stname=np.append(stname,"".join(np.array(auxstationname[i,:]).astype('str')))

				# Initial arrays
				# members, point, cycle/time, forecast_lead_time
				mwnd=np.zeros((size(members),size(stname),size(wlist),flt2.shape[0]),'f')*np.nan
				mhs=np.zeros((size(members),size(stname),size(wlist),flt2.shape[0]),'f')*np.nan
				mtm=np.zeros((size(members),size(stname),size(wlist),flt2.shape[0]),'f')*np.nan
				mdm=np.zeros((size(members),size(stname),size(wlist),flt2.shape[0]),'f')*np.nan
				# buoy data
				bwnd=np.zeros((size(stname),size(wlist),flt2.shape[0]),'f')*np.nan
				bhs=np.zeros((size(stname),size(wlist),flt2.shape[0]),'f')*np.nan
				btm=np.zeros((size(stname),size(wlist),flt2.shape[0]),'f')*np.nan
				bdm=np.zeros((size(stname),size(wlist),flt2.shape[0]),'f')*np.nan
				lat=np.zeros(size(stname),'f')*np.nan; lon=np.zeros(size(stname),'f')*np.nan

			else:
				if j==0:
					mtime = np.append(mtime,np.array([flt2[0]+sflt]).astype('double'),axis=0)

			awnd = np.array(ds2['wnd'].values[:,:]).T
			ahs = np.array(ds['hs'].values[:,:]).T
			atm = np.array(ds['tr'].values[:,:]).T
			adm = np.array(ds['th1m'].values[:,:]).T
			ds.close(); ds2.close(); del ds, ds2
			# Wind
			ignore,ind1,ind2 = np.intersect1d(mtime[t,:], flt2, return_indices=True)
			mwnd[j,:,t,:][:,ind1] = awnd[:,ind2]; del ind1, ind2
			# Wave
			ignore,ind1,ind2 = np.intersect1d(mtime[t,:], flt, return_indices=True)
			mhs[j,:,t,:][:,ind1] = ahs[:,ind2]
			mtm[j,:,t,:][:,ind1] = atm[:,ind2]
			mdm[j,:,t,:][:,ind1] = adm[:,ind2]
			del awnd,ahs,atm,adm,ind1,ind2,ignore

		print("completed time "+repr(t))

	print("completed ensemble member "+repr(j))

# help reading NDBC buoys, divided by year
yrange=np.array(np.arange(time.gmtime(mtime.min())[0],time.gmtime(mtime.min())[0]+1,1)).astype('int')
# Buoys
for b in range(0,size(stname)):

	ahs=[]
	try:
		awnd=[];ahs=[];atm=[];adm=[];atime=[]
		for y in yrange:
			f=nc.Dataset(ndbcp+"/"+stname[b]+"h"+repr(y)+".nc")
			awnd = np.append(awnd,f.variables['wind_spd'][:,0,0])
			ahs = np.append(ahs,f.variables['wave_height'][:,0,0])
			atm = np.append(atm,f.variables['average_wpd'][:,0,0])
			adm = np.append(adm,f.variables['mean_wave_dir'][:,0,0])
			atime = np.append(atime,np.array(f.variables['time'][:]).astype('double'))
			lat[b] = f.variables['latitude'][:]; lon[b] = f.variables['longitude'][:]
			f.close(); del f

	except:
		try:
			f=nc.Dataset(copernp+"/GL_TS_MO_"+stname[b]+".nc")
			awnd = np.nanmean(f.variables['WSPD'][:,:],axis=1)
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
						if np.any(awnd[indt[0]].mask==False):
							bwnd[b,t,i] = np.nanmean(awnd[indt[0]][awnd[indt[0]].mask==False])
						if np.any(atm[indt[0]].mask==False):
							btm[b,t,i] = np.nanmean(atm[indt[0]][atm[indt[0]].mask==False])
						if np.any(adm[indt[0]].mask==False):
							bdm[b,t,i] = np.nanmean(adm[indt[0]][adm[indt[0]].mask==False])

			print("counted "+repr(c)+" at "+stname[b])

	print("done "+stname[b])
	del ahs

# Simple quality-control (range)
ind=np.where((bwnd>80.)|(bwnd<0.0))
if size(ind)>0:
	bwnd[ind]=np.nan; del ind

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

ind=np.where((mwnd>80.)|(mwnd<0.0))
if size(ind)>0:
	mwnd[ind]=np.nan; del ind

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
ncfile = nc.Dataset("GEFSww3buoy_collocation_"+inidate+"to"+findate+".nc", "w", format=fnetcdf) 
ncfile.history="Collocation of GEFS WW3 point output ensemble forecast (table) and NDBC/Copernicus Buoys. Total of "+repr(bhs[bhs>0.].shape[0])+" observations or pairs model/observation."
# create  dimensions. 2 Dimensions
ncfile.createDimension('ensmembers', mhs.shape[0] )
ncfile.createDimension('buoypoints', bhs.shape[0] )
ncfile.createDimension('fcycletime', bhs.shape[1] )
ncfile.createDimension('forecastleadtime', bhs.shape[2] )
# create variables.
vt = ncfile.createVariable('time',np.dtype('float64').char,('fcycletime','forecastleadtime'))
vmembers = ncfile.createVariable('ensmembers',dtype('a25'),('ensmembers'))
vstname = ncfile.createVariable('buoyID',dtype('a25'),('buoypoints'))
vlat = ncfile.createVariable('latitude',np.dtype('float32').char,('buoypoints'))
vlon = ncfile.createVariable('longitude',np.dtype('float32').char,('buoypoints'))
#
vmwnd = ncfile.createVariable('model_wsp',np.dtype('float32').char,('ensmembers','buoypoints','fcycletime','forecastleadtime'))
vmhs = ncfile.createVariable('model_hs',np.dtype('float32').char,('ensmembers','buoypoints','fcycletime','forecastleadtime'))
vmtm = ncfile.createVariable('model_tm',np.dtype('float32').char,('ensmembers','buoypoints','fcycletime','forecastleadtime'))
vmdm = ncfile.createVariable('model_dm',np.dtype('float32').char,('ensmembers','buoypoints','fcycletime','forecastleadtime'))
vbwnd = ncfile.createVariable('buoy_wsp',np.dtype('float32').char,('buoypoints','fcycletime','forecastleadtime'))
vbhs = ncfile.createVariable('buoy_hs',np.dtype('float32').char,('buoypoints','fcycletime','forecastleadtime'))
vbtm = ncfile.createVariable('buoy_tm',np.dtype('float32').char,('buoypoints','fcycletime','forecastleadtime'))
vbdm = ncfile.createVariable('buoy_dm',np.dtype('float32').char,('buoypoints','fcycletime','forecastleadtime'))
# Assign units
vlat.units = 'degrees_north' ; vlon.units = 'degrees_east'
vt.units = 'seconds since 1970-01-01T00:00:00+00:00'
vmwnd.units='m/s'; vbwnd.units='m/s'
vmhs.units='m'; vbhs.units='m'
vmtm.units='s'; vbtm.units='s'
vmdm.units='degrees'; vbdm.units='degrees'
# Allocate Data
vt[:,:]=mtime[:,:]; vmembers[:]=members[:]; vstname[:]=stname[:]
vlat[:] = lat[:]; vlon[:] = lon[:]
vmwnd[:,:,:,:]=mwnd[:,:,:,:]
vmhs[:,:,:,:]=mhs[:,:,:,:]
vmtm[:,:,:,:]=mtm[:,:,:,:]
vmdm[:,:,:,:]=mdm[:,:,:,:]
vbwnd[:,:,:]=bwnd[:,:,:]
vbhs[:,:,:]=bhs[:,:,:]
vbtm[:,:,:]=btm[:,:,:]
vbdm[:,:,:]=bdm[:,:,:]
#
ncfile.close()
print(' ')
print('netcdf ok ')
