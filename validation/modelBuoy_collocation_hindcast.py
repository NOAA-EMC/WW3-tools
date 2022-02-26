# Collocation/pairing ww3 point output results with wave buoys.
# This code is restricted to ww3 hindcasts, where multiple files can be appended, creating a continuous time array.

import numpy as np
from matplotlib.mlab import *
from pylab import *
import xarray as xr
import netCDF4 as nc
import time
from time import strptime
from calendar import timegm
from pylab import *
import sys
import warnings; warnings.filterwarnings("ignore")
# netcdf format
fnetcdf="NETCDF4"

# Users have to check the buoys' names at the "f=nc.Dataset" lines below.

# Paths
# WW3 Model
mpath="/ww3runs/c00"
# NDBC buoys
ndbcp="/data/buoys/NDBC/wparam"
# Copernicus buoys
copernp="/data/buoys/Copernicus/wtimeseries"
# read list of ww3 files to be included in the collocation
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
		# list of station/buoy names
		if t==0:
			auxstationname=ds['station_name'].values[:,:]; stname=[]
			for i in range(0,auxstationname.shape[0]):
				stname=np.append(stname,"".join(np.array(auxstationname[i,:]).astype('str')))

		ahs = np.array(ds['hs'].values[:,:]).T
		atm = np.array(ds['tr'].values[:,:]).T
		adm = np.array(ds['th1m'].values[:,:]).T
		at = np.array(f.variables['time'][:]*24*3600 + timegm( strptime(np.str(f.variables['time'].units).split(' ')[2][0:4]+'01010000', '%Y%m%d%H%M') )).astype('double')
		ds.close();f.close(); del ds,f
		if t==0:
			mhs=np.copy(ahs)
			mtm=np.copy(atm)
			mdm=np.copy(adm)
			mtime=np.copy(at)
		else:
			mhs=np.append(mhs,ahs,axis=1)
			mtm=np.append(mtm,atm,axis=1)
			mdm=np.append(mdm,adm,axis=1)
			mtime=np.append(mtime,at)

		del ahs,atm,adm,at

# Buoys
bhs=np.zeros((size(stname),size(mtime)),'f')*np.nan
btm=np.zeros((size(stname),size(mtime)),'f')*np.nan
bdm=np.zeros((size(stname),size(mtime)),'f')*np.nan
lat=np.zeros(size(stname),'f')*np.nan; lon=np.zeros(size(stname),'f')*np.nan
for b in range(0,size(stname)):

	ahs=[]
	try:
		f=nc.Dataset(ndbcp+"/"+stname[b]+"h2016.nc")
		ahs = f.variables['wave_height'][:,0,0]
		atm = f.variables['average_wpd'][:,0,0]
		adm = f.variables['mean_wave_dir'][:,0,0]
		atime = np.array(f.variables['time'][:]).astype('double')
		lat[b] = f.variables['latitude'][:]; lon[b] = f.variables['longitude'][:]
		f.close(); del f
	except:
		try:
			f=nc.Dataset(copernp+"/GL_TS_MO_"+stname[b]+".nc")
			ahs = np.nanmean(f.variables['VHM0'][:,:],axis=1)
			atm = np.nanmean(f.variables['VTM02'][:,:],axis=1)
			adm = np.nanmean(f.variables['VMDR'][:,:],axis=1)
			atime = np.array(f.variables['TIME'][:]*24*3600 + timegm( strptime('195001010000', '%Y%m%d%H%M') )).astype('double')
			lat[b] = np.nanmean(f.variables['LATITUDE'][:]); lon[b] = np.nanmean(f.variables['LONGITUDE'][:])
			f.close(); del f
		except:
			ahs=[]

	else:
		if size(ahs>0):
			c=0
			for t in range(0,size(mtime)):
				indt=np.where(np.abs(atime-mtime[t])<1800.)
				if size(indt)>0:
					if np.any(ahs[indt[0]].mask==False):
						bhs[b,t] = np.nanmean(ahs[indt[0]][ahs[indt[0]].mask==False])
						c=c+1
					if np.any(atm[indt[0]].mask==False):
						btm[b,t] = np.nanmean(atm[indt[0]][atm[indt[0]].mask==False])
					if np.any(adm[indt[0]].mask==False):
						bdm[b,t] = np.nanmean(adm[indt[0]][adm[indt[0]].mask==False])

					del indt

			# print("counted "+repr(c)+" at "+stname[b])

	print("done "+stname[b])
	del ahs

# Save netcdf output file
lon[lon>180.]=lon[lon>180.]-360.
inidate=np.str(time.gmtime(mtime.min())[0])+np.str(time.gmtime(mtime.min())[1]).zfill(2)+np.str(time.gmtime(mtime.min())[2]).zfill(2)+np.str(time.gmtime(mtime.min())[3]).zfill(2)
findate=np.str(time.gmtime(mtime.max())[0])+np.str(time.gmtime(mtime.max())[1]).zfill(2)+np.str(time.gmtime(mtime.max())[2]).zfill(2)+np.str(time.gmtime(mtime.max())[3]).zfill(2)
ncfile = nc.Dataset("ww3buoy_collocation_"+inidate+"to"+findate+".nc", "w", format=fnetcdf) 
ncfile.history="Collocation of WW3 point output (table) and NDBC/Copernicus Buoys. Total of "+repr(bhs[bhs>0.].shape[0])+" observations or pairs model/observation."
# create  dimensions. 2 Dimensions
ncfile.createDimension('buoypoints', bhs.shape[0] )
ncfile.createDimension('time', bhs.shape[1] )
# create variables.
vt = ncfile.createVariable('time',np.dtype('float64').char,('time'))
vstname = ncfile.createVariable('buoyID',dtype('a25'),('buoypoints'))
vlat = ncfile.createVariable('latitude',np.dtype('float32').char,('buoypoints'))
vlon = ncfile.createVariable('longitude',np.dtype('float32').char,('buoypoints'))
#
vmhs = ncfile.createVariable('model_hs',np.dtype('float32').char,('buoypoints','time'))
vmtm = ncfile.createVariable('model_tm',np.dtype('float32').char,('buoypoints','time'))
vmdm = ncfile.createVariable('model_dm',np.dtype('float32').char,('buoypoints','time'))
vbhs = ncfile.createVariable('buoy_hs',np.dtype('float32').char,('buoypoints','time'))
vbtm = ncfile.createVariable('buoy_tm',np.dtype('float32').char,('buoypoints','time'))
vbdm = ncfile.createVariable('buoy_dm',np.dtype('float32').char,('buoypoints','time'))
# Assign units
vlat.units = 'degrees_north' ; vlon.units = 'degrees_east'
vt.units = 'seconds since 1970-01-01T00:00:00+00:00'
vmhs.units='m'; vbhs.units='m'
vmtm.units='s'; vbtm.units='s'
vmdm.units='degrees'; vbdm.units='degrees'
# Allocate Data
vt[:]=mtime[:]; vstname[:]=stname[:]
vlat[:] = lat[:]; vlon[:] = lon[:]
vmhs[:,:]=mhs[:,:]
vmtm[:,:]=mtm[:,:]
vmdm[:,:]=mdm[:,:]
vbhs[:,:]=bhs[:,:]
vbtm[:,:]=btm[:,:]
vbdm[:,:]=bdm[:,:]
#
ncfile.close()
print(' ')
print('netcdf ok ')

