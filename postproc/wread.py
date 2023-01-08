#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
wread.py

VERSION AND LAST UPDATE:
 v1.0  04/04/2022
 v1.1  01/06/2023

PURPOSE:
 Group of python functions to Read Wave data: 
  WAVEWATCHIII results, and NDBC and Copernicus buoys.
 Prefix meaning:
 tseriesnc = time series (table of integrated parameters versus time).
 spec = wave spectrum.
 Users can import as a standard python function, and use it accordingly:
 For example:
  import wread
  wread.tseriesnc_ww3(filename.nc,stationID)
 Users can help() each function to obtain information about inputs/outputs
  help(wread.tseriesnc_ww3)

USAGE:
 functions
   mask
   cyclonemap
   tseriesnc_ndbc
   tseriesnc_copernicus
   tseriesnc_ww3
   bull
   bull_tar
   ts
   station_tar
   spec_ndbc
   spec_ww3
 Explanation for each function is contained in the headers

OUTPUT:
 numpy arrays. Description of variables is contained in the header 
  of each function 

DEPENDENCIES:
 See dependencies.py and the imports below.

AUTHOR and DATE:
 04/04/2022: Ricardo M. Campos, first version.
 01/06/2023: Ricardo M. Campos, new file formats added.

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import matplotlib
matplotlib.use('Agg') # for backend plots, not for rendering in a window
import time
from time import strptime
from calendar import timegm
import xarray as xr
import netCDF4 as nc
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import sys
import pandas as pd
from matplotlib import ticker
# import pickle
import sys
import warnings; warnings.filterwarnings("ignore")


# FIELDS

def mask(*args):
	'''
	Read gridmask netcdf file generated with prepGridMask.py
	Input: file name (example: gridInfo_GEFSv12.nc)
	Output: dictionary containing the arrays and string names
	'''
	if len(args) == 1:
		fname=np.str(args[0])
	else:
		sys.exit(' Too many inputs')

	print("  reading ww3_tools mask ...")
	try:
		f=nc.Dataset(fname)
		# build dictionary
		result={'latitude':np.array(f.variables['latitude'][:]),'longitude':np.array(f.variables['longitude'][:]),'mask':np.array(f.variables['mask'][:,:])}
	except:
		sys.exit(" Cannot open "+fname)
	else:
		if 'distcoast' in f.variables.keys():
			result['distcoast'] = np.array(f.variables['distcoast'][:,:])
		if 'depth' in f.variables.keys():
			result['depth'] = np.array(f.variables['depth'][:,:])	
		if 'GlobalOceansSeas' in f.variables.keys():
			result['GlobalOceansSeas'] = np.array(f.variables['GlobalOceansSeas'][:,:])		
		if 'HighSeasMarineZones' in f.variables.keys():
			result['HighSeasMarineZones'] = np.array(f.variables['HighSeasMarineZones'][:,:])	
		if 'names_GlobalOceansSeas' in f.variables.keys():
			result['names_GlobalOceansSeas'] = f.variables['names_GlobalOceansSeas'][:]
		if 'names_HighSeasMarineZones' in f.variables.keys():
			result['names_HighSeasMarineZones'] = f.variables['names_HighSeasMarineZones'][:]

		f.close(); del f
		print("  GridMask Read. "+fname)
		return result
		del result

def cyclonemap(*args):
	'''
	Read cyclonemap netcdf file generated with procyclmap.py
	Input: file name (example: CycloneMap2020.nc)
	Output: dictionary containing the arrays and string names
	'''
	if len(args) == 1:
		fname=np.str(args[0])
	else:
		sys.exit(' Too many inputs')

	print("  reading ww3_tools cyclonemap ...")
	try:
		f=nc.MFDataset(fname, aggdim='time')
		at=f.variables['time'][:]; adate=[]
		for j in range(0,at.shape[0]):
			adate=np.append(adate,date2num(datetime.datetime(time.gmtime(at[j])[0],time.gmtime(at[j])[1],time.gmtime(at[j])[2],time.gmtime(at[j])[3],time.gmtime(at[j])[4])))		
		# --------
		# build dictionary
		result={'latitude':np.array(f.variables['lat'][:]),'longitude':np.array(f.variables['lon'][:]),
		'time':np.array(at).astype('double'),'date':np.array(adate).astype('double'),
		'cmap':f.variables['cmap'], 'info':np.str(f.info), 'netcdf':f}
		# it does not allocate the data of cmap using [:,:,:] yet as it can take a lot of data/memory and time.
	except:
		sys.exit(" Cannot open "+fname)
	else:
		print("    CycloneInfo Read. "+fname)
		return result
		del result


# TIME-SERIES

# Observations NDBC, netcdf format
def tseriesnc_ndbc(*args):
	'''
	Observations NDBC, time series/table, netcdf format
	Input: file name (example: 46047h2016.nc)
	Output: dictionary containing the arrays: time(seconds since 1970),time(datetime64),lat,lon, 
	  and arrays sst,mslp,dwp,tmp,gst,wsp,wdir,hs,tm,tp,dm
	'''
	if len(args) == 1:
		fname=np.str(args[0])
	elif len(args) > 1:
		sys.exit(' Too many inputs')

	try:
		ds = xr.open_dataset(fname); f=nc.Dataset(fname)
	except:
		sys.exit(" Cannot open "+fname)
	else:
		btm = f.variables['average_wpd'][:,0,0]; btp = f.variables['dominant_wpd'][:,0,0]
		btime = np.array(f.variables['time'][:]).astype('double')
		f.close(); del f
		bsst = ds['sea_surface_temperature'].values[:,0,0]
		bmslp = ds['air_pressure'].values[:,0,0]
		bdwp = ds['dewpt_temperature'].values[:,0,0]
		btmp = ds['air_temperature'].values[:,0,0]
		bgst = ds['gust'].values[:,0,0]	

		if 'wind_spd' in ds.keys():
			bwsp = ds['wind_spd'].values[:,0,0]
			try:
				from urllib.request import urlopen
				url = "https://www.ndbc.noaa.gov/station_page.php?station="+fname.split('h')[0]
				page = urlopen(url)
				html_bytes = page.read()
				html = html_bytes.decode("utf-8")
			except:
				anh=4.0 # assuming most of anemometer heights are between 3.7 to 4.1.
				print('Information about the Anemometer height, for wind speed conversion to 10m, could not be obtained.')
			else:
				if "Anemometer height" in html:
					anh=np.float(html.split('Anemometer height')[1][0:15].split(':</b>')[1].split('m')[0])
				else:
					print('Information about the Anemometer height, for wind speed conversion to 10m, could not be found.')
					anh=4.0 # assuming most of anemometer heights are between 3.7 to 4.1.

				del url,page,html_bytes,html

			# convert wind speed to 10 meters (DNVGL C-205, confirmed by https://onlinelibrary.wiley.com/doi/pdf/10.1002/er.6382)
			bwsp =  np.copy(((10./anh)**(0.1)) * bwsp)

		bwdir = ds['wind_dir'].values[:]
		bhs = ds['wave_height'].values[:,0,0]
		bdm = ds['mean_wave_dir'].values[:,0,0]

		# Automatic and basic Quality Control
		bsst[np.abs(bsst)>70]=np.nan
		bmslp[(bmslp<500)|(bmslp>1500)]=np.nan
		bdwp[np.abs(bdwp)>80]=np.nan
		btmp[np.abs(btmp)>80]=np.nan
		bgst[(bgst<0)|(bgst>200)]=np.nan
		bwsp[(bwsp<0)|(bwsp>150)]=np.nan
		bwdir[(bwdir<-180)|(bwdir>360)]=np.nan
		bhs[(bhs<0)|(bhs>30)]=np.nan
		btm[(btm<0)|(btm>40)]=np.nan
		btp[(btp<0)|(btp>40)]=np.nan
		bdm[(bdm<-180)|(bdm>360)]=np.nan

		result={'latitude':np.array(ds['latitude'].values[:]),'longitude':np.array(ds['longitude'].values[:]),
		'time':btime,'date':ds['time'].values[:],
		'sst':bsst, 'mslp':bmslp, 'dewpt_temp':bdwp,
		'air_temp':btmp, 'gust':bgst, 'wind_spd':bwsp,	
		'wind_dir':bwdir, 'hs':bhs, 'tm':btm,	
		'tp':btp, 'dm':bdm, 'tm':btm}

		return result
		ds.close()
		del ds,btime,bsst,bmslp,bdwp,btmp,bgst,bwsp,bwdir,bhs,btm,btp,bdm

# Observations Copernicus, netcdf format
def tseriesnc_copernicus(*args):
	'''
	Observations NDBC, time series/table, netcdf format
	Input: file name (example: 46047h2016.nc)
	Output: dictionary containing the arrays: time(seconds since 1970),time(datetime64),lat,lon, 
	  and arrays with the environmental variables available.
	'''
	if len(args) == 1:
		fname=np.str(args[0])
	elif len(args) > 1:
		sys.exit(' Too many inputs')

	try:
		ds = xr.open_dataset(fname); f=nc.Dataset(fname)
	except:
		sys.exit(" Cannot open "+fname)
	else:
		btime = np.array(f.variables['TIME'][:]*24*3600 + timegm( strptime('195001010000', '%Y%m%d%H%M') )).astype('double')
		f.close(); del f
		blat = np.nanmean(ds['LATITUDE'].values[:])
		blon = np.nanmean(ds['LONGITUDE'].values[:])
		# dictionary
		result={'latitude':np.array(blat),'longitude':np.array(blon),
		'time':btime,'date':ds['TIME'].values[:]}

		if 'DEPH' in ds.keys():		
			bdepth = np.nanmean(ds['DEPH'].values[:,:],axis=1) # Depth 
			result['depth']=np.array(bdepth)

		if 'VHM0' in ds.keys():
			bhs = np.nanmean(ds['VHM0'].values[:,:],axis=1) # Hs
			bhs[(bhs<0)|(bhs>30)]=np.nan
			result['hs']=np.array(bhs)
		elif 'VGHS' in ds.keys():
			bhs = np.nanmean(ds['VGHS'].values[:,:],axis=1) # Hs
			bhs[(bhs<0)|(bhs>30)]=np.nan
			result['hs']=np.array(bhs)

		if 'VAVH' in ds.keys():  
			bvavh = np.nanmean(ds['VAVH'].values[:,:],axis=1) # H 1/3 vavh
			bvavh[(bvavh<0)|(bvavh>30)]=np.nan
			result['hs_vavh']=np.array(bvavh)

		if 'VZMX' in ds.keys():
			bhmax = np.nanmean(ds['VZMX'].values[:,:],axis=1) # Hmax
			bhmax[(bhmax<0)|(bhmax>40)]=np.nan
			result['hmax']=np.array(bhmax)

		if 'VTM02' in ds.keys():
			btm = np.nanmean(ds['VTM02'].values[:,:],axis=1) # Tm
			btm[(btm<0)|(btm>40)]=np.nan
			result['tm']=np.array(btm)
		elif 'VGTA' in ds.keys():
			btm = np.nanmean(ds['VGTA'].values[:,:],axis=1) # Tm
			btm[(btm<0)|(btm>40)]=np.nan
			result['tm']=np.array(btm)

		if 'VTPK' in ds.keys():
			btp = np.nanmean(ds['VTPK'].values[:,:],axis=1) # Tp
			btp[(btp<0)|(btp>40)]=np.nan
			result['tp']=np.array(btp)

		if 'TEMP' in ds.keys():		
			bsst = np.nanmean(ds['TEMP'].values[:,:],axis=1) # SST 
			bsst[np.abs(bsst)>70]=np.nan
			result['sst']=np.array(bsst)

		if 'ATMS' in ds.keys():
			bmslp = np.nanmean(ds['ATMS'].values[:,:],axis=1) # Pressure
			bmslp[(bmslp<500)|(bmslp>1500)]=np.nan
			result['mslp']=np.array(bmslp)

		if 'DEWT' in ds.keys():
			bdwp = np.nanmean(ds['DEWT'].values[:,:],axis=1) # Dewpoint
			bdwp[np.abs(bdwp)>80]=np.nan
			result['dewpt_temp']=np.array(bdwp)

		if 'DRYT' in ds.keys():
			btmp = np.nanmean(ds['DRYT'].values[:,:],axis=1) # air temperature
			btmp[np.abs(btmp)>80]=np.nan
			result['air_temp']=np.array(btmp)

		if 'GSPD' in ds.keys():
			bgst = np.nanmean(ds['GSPD'].values[:,:],axis=1) # gust
			bgst=np.copy(((10./4.0)**(0.1))*bgst) # conversion to 10m, approximation
			bgst[(bgst<0)|(bgst>200)]=np.nan
			result['gust']=np.array(bgst)

		if 'WSPD' in ds.keys():
			bwsp = np.nanmean(ds['WSPD'].values[:,:],axis=1) # wind speed		
			bwsp=np.copy(((10./4.0)**(0.1))*bwsp) # conversion to 10m, approximation
			bwsp[(bwsp<0)|(bwsp>150)]=np.nan
			result['wind_spd']=np.array(bwsp)

		if 'WDIR' in ds.keys():
			bwdir = np.nanmean(ds['WDIR'].values[:,:],axis=1) # wind direction		
			bwdir[(bwdir<-180)|(bwdir>360)]=np.nan
			result['wind_dir']=np.array(bwdir)

		if 'VCMX' in ds.keys():
			bhcmax = np.nanmean(ds['VCMX'].values[:,:],axis=1) # Hcrest max
			bhcmax[(bhcmax<0)|(bhcmax>40)]=np.nan
			result['hc_max']=np.array(bhcmax)

		if 'VMDR' in ds.keys():	
			bdm = np.nanmean(ds['VMDR'].values[:,:],axis=1) # Mean direction
			bdm[(bdm<-180)|(bdm>360)]=np.nan
			result['dm']=np.array(bdm)

		if 'VPED' in ds.keys():
			bdp = np.nanmean(ds['VPED'].values[:,:],axis=1) # Peak direction
			bdp[(bdp<-180)|(bdp>360)]=np.nan
			result['dp']=np.array(bdp)

		return result
		ds.close(); del ds


# WAVEWATCH III point output, netcdf format
def tseriesnc_ww3(*args):
	'''
	WAVEWATCH III, time series/table, netcdf format
	Input:  file name (example: ww3gefs.20160928_tab.nc), and station name (example: 41002)
	Output: dictionary containing the arrays: time(seconds since 1970),time(datetime64),lat,lon, 
	  and arrays with the wave variables available.
	'''
	if len(args) == 2:
		fname=np.str(args[0]); stname=np.str(args[1])
	elif len(args) < 2 :
		sys.exit(' Two inputs are required: file name and station name')
	elif len(args) > 2:
		sys.exit(' Too many inputs')

	try:
		ds = xr.open_dataset(fname); f=nc.Dataset(fname)
	except:
		sys.exit(" Cannot open "+fname)
	else:
		mtime = np.array(f.variables['time'][:]*24*3600 + timegm( strptime(np.str(f.variables['time'].units).split(' ')[2][0:4]+'01010000', '%Y%m%d%H%M') )).astype('double')
		f.close(); del f
	
		auxstationname=ds['station_name'].values[:,:]; stationname=[]
		for i in range(0,auxstationname.shape[0]):
			stationname=np.append(stationname,"".join(np.array(auxstationname[i,:]).astype('str')))

		inds=np.where(stationname[:]==stname)
		if size(inds)>0:
			inds=np.int(inds[0][0]); stname=np.str(stationname[inds])
		else:
			sys.exit(' Station '+stname+' not included in the ww3 output file, or wrong station ID')

		mlat = np.nanmean(ds['latitude'].values[:,inds])
		mlon = np.nanmean(ds['longitude'].values[:,inds])
		# dictionary
		result={'latitude':np.array(mlat),'longitude':np.array(mlon),
		'time':mtime,'date':ds['time'].values[:]}	

		if 'hs' in ds.keys():
			mhs = ds['hs'].values[:,inds]
			mhs[(mhs<0)|(mhs>30)]=np.nan
			result['hs']=np.array(mhs)
		if 'fp' in ds.keys():
			mtp = np.zeros(mhs.shape[0],'f')*np.nan
			mtp[ds['fp'].values[:,inds]>0.0] = 1./ds['fp'].values[:,inds][ds['fp'].values[:,inds]>0.0]
			mtp[(mtp<0)|(mtp>40)]=np.nan
			result['tp']=np.array(mtp)
		if 'tr' in ds.keys():
			mtm = ds['tr'].values[:,inds]
			mtm[(mtm<0)|(mtm>40)]=np.nan
			result['tm']=np.array(mtm)
		if 'th1p' in ds.keys():
			mdp = ds['th1p'].values[:,inds]
			mdp[(mdp<-180)|(mdp>360)]=np.nan
			result['dp']=np.array(mdp)
		if 'th1m' in ds.keys():
			mdm = ds['th1m'].values[:,inds]
			mdm[(mdm<-180)|(mdm>360)]=np.nan
			result['dm']=np.array(mdm)
		if 'sth1m' in ds.keys():
			result['spr']=np.array(ds['sth1m'].values[:,inds])
		if 'lm' in ds.keys():
			result['lm']=np.array(ds['lm'].values[:,inds])
		if 'sth1p' in ds.keys():
			result['spr_dp']=np.array(ds['sth1p'].values[:,inds])

		return result
		ds.close(); del ds


# Operational WW3 formats

def bull(*args):
	'''
	WAVEWATCH III, bull operational point output, see https://www.ftp.ncep.noaa.gov/data/nccf/com/
	Input:  file name (example: gefs.wave.41004.bull)
	Output: dictionary containing:
	  time(seconds since 1970),time(datetime64),lat,lon,station name; Arrays: hs, tp, and dp (gfs only)
	'''
	if len(args) == 1:
		fname=np.str(args[0])
	elif len(args) > 1:
		sys.exit(' Too many inputs')

	# confirm format
	if np.str(fname).split('/')[-1].split('.')[-1]=='bull':
		print("  reading ww3 bull file ...")
		at=[]; adate=[]; ahs=[]; atp=[]; adp=[]
		stname=np.str(fname).split('/')[-1].split('.')[-2]

		try:
			tfile = open(fname, 'r'); lines = tfile.readlines()
		except:
			sys.exit('   Cannot open '+fname)
		else:

			if 'gfs' in np.str(fname).split('/')[-1]:
				iauxhs=[24,30];iauxtp=[30,34];iauxdp=[35,38]

				# lat / lon
				auxpos=np.str(lines[0]).replace("b'","").split('(')[1]
				if auxpos[5]=='N':
					alat=np.float(auxpos[0:5])
				else:
					alat=-1.*np.float(auxpos[0:5])

				if auxpos[13]=='E':
					alon=np.float(auxpos[7:13])
				else:
					alon=-1.*np.float(auxpos[7:13])

				# time ----
				auxdate = np.str(lines[2]).split(':')[1].split('UTC')[0][1::]
				auxt = np.double(timegm( strptime(  auxdate[0:8]+' '+auxdate[9:11]+'00', '%Y%m%d %H%M') ))
				year = np.int(time.gmtime(auxt)[0]); month = np.int(time.gmtime(auxt)[1])
				pday=0
				for j in range(7,size(lines)-8):
					day=np.int(lines[j][3:5]); hour=np.int(lines[j][6:8])
					if day<pday:
						if month<12:
							month=month+1
						else:
							month=1; year=year+1

					at=np.append(at,np.double(timegm( strptime( repr(year)+np.str(month).zfill(2)+np.str(day).zfill(2)+' '+np.str(hour).zfill(2)+'00', '%Y%m%d %H%M') )))
					pday=np.copy(day)

				del hour,day,month,year
				for j in range(0,at.shape[0]):
					adate=np.append(adate,date2num(datetime.datetime(time.gmtime(at[j])[0],time.gmtime(at[j])[1],time.gmtime(at[j])[2],time.gmtime(at[j])[3],time.gmtime(at[j])[4])))		
				# --------

				for j in range(7,size(lines)-8):
					if len(lines[j][iauxhs[0]:iauxhs[1]].replace(' ',''))>0:
						ahs=np.append(ahs,np.float(lines[j][10:15]))
						auxhs=np.array([np.float(lines[j][iauxhs[0]:iauxhs[1]])])
						for k in range(1,4):		
							if len(np.str(lines[j][int(iauxhs[0]+18*k):int(iauxhs[1]+18*k)]).replace(' ', '')):
								auxhs=np.append(auxhs,np.float(lines[j][int(iauxhs[0]+18*k):int(iauxhs[1]+18*k)]))

						auxtp=np.array([np.float(lines[j][iauxtp[0]:iauxtp[1]])])
						for k in range(1,4):		
							if len(np.str(lines[j][int(iauxtp[0]+18*k):int(iauxtp[1]+18*k)]).replace(' ', '')):
								auxtp=np.append(auxtp,np.float(lines[j][int(iauxtp[0]+18*k):int(iauxtp[1]+18*k)]))

						auxdp=np.array([np.float(lines[j][iauxdp[0]:iauxdp[1]])])
						for k in range(1,4):		
							if len(np.str(lines[j][int(iauxdp[0]+18*k):int(iauxdp[1]+18*k)]).replace(' ', '')):
								auxdp=np.append(auxdp,np.float(lines[j][int(iauxdp[0]+18*k):int(iauxdp[1]+18*k)]))

						indaux=np.nanmin(np.where(auxhs==np.nanmax(auxhs))[0])
						atp=np.append(atp,np.float(auxtp[indaux]))
						adp=np.append(adp,np.float(auxdp[indaux]))
						del indaux,auxhs,auxtp,auxdp
					else:
						ahs=np.append(ahs,np.nan)
						atp=np.append(atp,np.nan)
						adp=np.append(adp,np.nan)

				# build dictionary
				result={'latitude':alat,'longitude':alon,'station_name':stname,
				'time':np.array(at).astype('double'),'date':np.array(adate).astype('double'),
				'hs':np.array(ahs),'tp':np.array(atp),'dp':np.array(adp)}

				del adp

			elif 'gefs' in np.str(fname).split('/')[-1]:
				iauxhs=[10,15];iauxtp=[28,33]

				# lat / lon
				auxpos=np.str(lines[1]).split('(')[1].split('N')
				alat=np.float(auxpos[0])
				alon=np.float(auxpos[1].split('W')[0])

				# time ----
				auxdate = np.str(lines[3]).split(':')[1].split('UTC')[0][1::]
				auxt = np.double(timegm( strptime(  auxdate[0:8]+' '+auxdate[10:12]+'00', '%Y%m%d %H%M') ))
				year = np.int(time.gmtime(auxt)[0]); month = np.int(time.gmtime(auxt)[1])
				pday=0
				for j in range(9,size(lines)-8):
					day=np.int(lines[j][2:4]); hour=np.int(lines[j][5:7])
					if day<pday:
						if month<12:
							month=month+1
						else:
							month=1; year=year+1

					at=np.append(at,np.double(timegm( strptime( repr(year)+np.str(month).zfill(2)+np.str(day).zfill(2)+' '+np.str(hour).zfill(2)+'00', '%Y%m%d %H%M') )))
					pday=np.copy(day)

				del hour,day,month,year
				for j in range(0,at.shape[0]):
					adate=np.append(adate,date2num(datetime.datetime(time.gmtime(at[j])[0],time.gmtime(at[j])[1],time.gmtime(at[j])[2],time.gmtime(at[j])[3],time.gmtime(at[j])[4])))		
				# --------

				ahs=[]; atp=[]
				for j in range(9,size(lines)-8):
					if len(lines[j][iauxhs[0]:iauxhs[1]].replace(' ',''))>0:
						ahs=np.append(ahs,np.float(lines[j][iauxhs[0]:iauxhs[1]]))
						atp=np.append(atp,np.float(lines[j][iauxtp[0]:iauxtp[1]]))
					else:
						ahs=np.append(ahs,np.nan)
						atp=np.append(atp,np.nan)

				# build dictionary
				result={'time':np.array(at).astype('double'),'date':np.array(adate).astype('double'),
				'latitude':alat,'longitude':alon,'station_name':stname,				
				'hs':np.array(ahs),'tp':np.array(atp)}

			print("   Model data read, "+fname+", bull format.")
			return result
			del result,alat,alon,at,adate,ahs,atp,tfile,lines
	else:
		sys.exit(" Skipped file "+fname+" Not bull_tar format.")


def bull_tar(*args):
	'''
	WAVEWATCH III, bull_tar operational point output, see https://www.ftp.ncep.noaa.gov/data/nccf/com/
	Input:  file name (example: gfswave.t00z.bull_tar)
	Output: dictionary containing:
	  time(seconds since 1970),time(datetime64),lat,lon,station names; Arrays: hs, tp, and dp (gfs only)
	'''
	if len(args) == 1:
		fname=np.str(args[0])
	elif len(args) > 1:
		sys.exit(' Too many inputs')

	# confirm file format
	if fname.split('/')[-1].split('.')[-1]=='bull_tar':
		print("  reading ww3 bull_tar file ...")
		import tarfile
		stname=[]

		try:
			tar = tarfile.open(fname)
		except:
			sys.exit('   Cannot open '+fname)
		else:
			at=[]; adate=[]; alat=[]; alon=[]

			if 'gfs' in np.str(fname).split('/')[-1]:
				iauxhs=[24,30];iauxtp=[30,34];iauxdp=[35,38]

				for t in range(0,size(tar.getmembers())):
					# station names
					stname=np.append(stname,np.str(np.str(tar.getmembers()[t].name).split('/')[-1]).split('/')[-1].split('.')[-2])

					try:
						tfile=tar.extractfile(tar.getmembers()[t]); lines = tfile.readlines()
					except:
						print("   Cannot open "+tar.getmembers()[t].name)
					else:

						# lat / lon
						auxpos=np.str(lines[0]).replace("b'","").split('(')[1]
						if auxpos[5]=='N':
							alat=np.append(alat,np.float(auxpos[0:5]))
						else:
							alat=np.append(alat,-1.*np.float(auxpos[0:5]))

						if auxpos[13]=='E':
							alon=np.append(alon,np.float(auxpos[7:13]))
						else:
							alon=np.append(alon,-1.*np.float(auxpos[7:13]))

						if t==0:				
							# time array ----
							auxdate = np.str(lines[2]).split(':')[1].split('UTC')[0][1::]
							auxt = np.double(timegm( strptime(  auxdate[0:8]+' '+auxdate[9:11]+'00', '%Y%m%d %H%M') ))
							year = np.int(time.gmtime(auxt)[0]); month = np.int(time.gmtime(auxt)[1])
							pday=0
							for j in range(7,size(lines)-8):
								auxlines = np.str(lines[j]).replace("b'","")
								day=np.int(auxlines[3:5]); hour=np.int(auxlines[6:8]); del auxlines
								if day<pday:
									if month<12:
										month=month+1
									else:
										month=1; year=year+1

								at=np.append(at,np.double(timegm( strptime( repr(year)+np.str(month).zfill(2)+np.str(day).zfill(2)+' '+np.str(hour).zfill(2)+'00', '%Y%m%d %H%M') )))
								pday=np.copy(day)

							del hour,day,month,year
							for j in range(0,at.shape[0]):
								adate=np.append(adate,date2num(datetime.datetime(time.gmtime(at[j])[0],time.gmtime(at[j])[1],time.gmtime(at[j])[2],time.gmtime(at[j])[3],time.gmtime(at[j])[4])))

							# --------
							ahs=np.zeros((size(tar.getmembers()),at.shape[0]),'f')*np.nan
							atp=np.zeros((size(tar.getmembers()),at.shape[0]),'f')*np.nan
							adp=np.zeros((size(tar.getmembers()),at.shape[0]),'f')*np.nan

						auxhs=[]; auxtp=[]; auxdp=[]
						for j in range(7,size(lines)-8):
							auxlines = np.str(lines[j]).replace("b'","")
							if len(auxlines[iauxhs[0]:iauxhs[1]].replace(' ',''))>0:
								auxhs=np.append(auxhs,np.float(auxlines[10:15]))
								fuxhs=np.array([np.float(auxlines[iauxhs[0]:iauxhs[1]])])
								for k in range(1,4):		
									if len(np.str(auxlines[int(iauxhs[0]+18*k):int(iauxhs[1]+18*k)]).replace(' ', '')):
										fuxhs=np.append(fuxhs,np.float(auxlines[int(iauxhs[0]+18*k):int(iauxhs[1]+18*k)]))

								fuxtp=np.array([np.float(auxlines[iauxtp[0]:iauxtp[1]])])
								for k in range(1,4):		
									if len(np.str(auxlines[int(iauxtp[0]+18*k):int(iauxtp[1]+18*k)]).replace(' ', '')):
										fuxtp=np.append(fuxtp,np.float(auxlines[int(iauxtp[0]+18*k):int(iauxtp[1]+18*k)]))

								fuxdp=np.array([np.float(auxlines[iauxdp[0]:iauxdp[1]])])
								for k in range(1,4):		
									if len(np.str(auxlines[int(iauxdp[0]+18*k):int(iauxdp[1]+18*k)]).replace(' ', '')):
										fuxdp=np.append(fuxdp,np.float(auxlines[int(iauxdp[0]+18*k):int(iauxdp[1]+18*k)]))

								indaux=np.nanmin(np.where(fuxhs==np.nanmax(fuxhs))[0])
								auxtp=np.append(auxtp,np.float(fuxtp[indaux]))
								auxdp=np.append(auxdp,np.float(fuxdp[indaux]))
								del indaux,fuxhs,fuxtp,fuxdp
							else:
								auxhs=np.append(auxhs,np.nan)
								auxtp=np.append(auxtp,np.nan)
								auxdp=np.append(auxdp,np.nan)

						if ahs.shape[1]==auxhs.shape[0]:
							ahs[t,:]=np.array(auxhs)
							atp[t,:]=np.array(auxtp)
							adp[t,:]=np.array(auxdp)
						else:
							print("   Time duration of "+tar.getmembers()[t]+" (in "+fname+") do not match the other stations. Mantained NaN.")

						del auxhs,auxtp,auxdp,tfile,lines

				# build dictionary			
				result={'time':np.array(at).astype('double'),'date':np.array(adate).astype('double'),
				'latitude':np.array(alat),'longitude':np.array(alon),'station_name':np.array(stname),
				'hs':np.array(ahs),'tp':np.array(atp),'dp':np.array(adp)}

				del adp

			elif 'gefs' in np.str(fname).split('/')[-1]:
				iauxhs=[10,15];iauxtp=[28,33]
				
				for t in range(0,size(tar.getmembers())):
					# station names
					stname=np.append(stname,np.str(np.str(tar.getmembers()[t].name).split('/')[-1]).split('/')[-1].split('.')[-2])

					try:
						tfile=tar.extractfile(tar.getmembers()[t]); lines = tfile.readlines()
					except:
						print("   Cannot open "+tar.getmembers()[t].name)
					else:
						# lat / lon
						auxpos=np.str(lines[1]).replace("b'","").split('(')[1].split('N')
						alat=np.append(alat,np.float(auxpos[0]))
						alon=np.append(alon,np.float(auxpos[1].split('W')[0]))

						if t==0:
							# time array ----
							auxdate = np.str(lines[3]).split(':')[1].split('UTC')[0][1::]
							auxt = np.double(timegm( strptime(  auxdate[0:8]+' '+auxdate[10:12]+'00', '%Y%m%d %H%M') ))
							year = np.int(time.gmtime(auxt)[0]); month = np.int(time.gmtime(auxt)[1])
							pday=0
							for j in range(9,size(lines)-8):
								auxlines = np.str(lines[j]).replace("b'","")
								day=np.int(auxlines[2:4]); hour=np.int(auxlines[5:7]); del auxlines
								if day<pday:
									if month<12:
										month=month+1
									else:
										month=1; year=year+1

								at=np.append(at,np.double(timegm( strptime( repr(year)+np.str(month).zfill(2)+np.str(day).zfill(2)+' '+np.str(hour).zfill(2)+'00', '%Y%m%d %H%M') )))
								pday=np.copy(day)

							del hour,day,month,year
							for j in range(0,at.shape[0]):
								adate=np.append(adate,date2num(datetime.datetime(time.gmtime(at[j])[0],time.gmtime(at[j])[1],time.gmtime(at[j])[2],time.gmtime(at[j])[3],time.gmtime(at[j])[4])))

							# --------
							ahs=np.zeros((size(tar.getmembers()),at.shape[0]),'f')*np.nan
							atp=np.zeros((size(tar.getmembers()),at.shape[0]),'f')*np.nan

						auxhs=[]; auxtp=[]
						for j in range(9,size(lines)-8):
							auxlines = np.str(lines[j]).replace("b'","")
							if len(auxlines[iauxhs[0]:iauxhs[1]].replace(' ',''))>0:
								auxhs=np.append(auxhs,np.float(auxlines[iauxhs[0]:iauxhs[1]]))
								auxtp=np.append(auxtp,np.float(auxlines[iauxtp[0]:iauxtp[1]]))
							else:
								auxhs=np.append(auxhs,np.nan)
								auxtp=np.append(auxtp,np.nan)

							del auxlines

						if ahs.shape[1]==auxhs.shape[0]:
							ahs[t,:]=np.array(auxhs)
							atp[t,:]=np.array(auxtp)
						else:
							print("   Time duration of "+tar.getmembers()[t]+" (in "+fname+") do not match the other stations. Mantained NaN.")
	
						del auxhs,auxtp,tfile,lines

				# build dictionary			
				result={'time':np.array(at).astype('double'),'date':np.array(adate).astype('double'),
				'latitude':np.array(alat),'longitude':np.array(alon),'station_name':np.array(stname),
				'hs':np.array(ahs),'tp':np.array(atp)}

			print("   Model data read, "+fname+", bull_tar format.")
			return result
			del result,tar,alat,alon,ahs,atp,at,adate
	else:
		sys.exit(" Skipped file "+fname+" Not bull_tar format.")


def ts(*args):
	'''
	WAVEWATCH III, ts operational point output, see https://www.ftp.ncep.noaa.gov/data/nccf/com/
	Input:  file name (example: gefs.wave.41004.ts)
	Output: dictionary containing:
	  time(seconds since 1970),time(datetime64),station name; Arrays: hs, hs_spr, tp (glwu or gefs)
	'''
	if len(args) == 1:
		fname=np.str(args[0])
	elif len(args) > 1:
		sys.exit(' Too many inputs')

	# confirm format
	if np.str(fname).split('/')[-1].split('.')[-1]=='ts':
		print("  reading ww3 ts file ...")
		stname=np.str(fname).split('/')[-1].split('.')[-2]
		try:
			tfile = pd.read_csv(fname,skiprows=2); lines = tfile.values[:,0]
		except:
			sys.exit('   Cannot open '+fname)
		else:

			if 'gefs' in np.str(fname).split('/')[-1]:
				# gefs lakes ww3 format
				at=[]; adate=[]; ahs=[]; ahspr=[]; atp=[]
				for j in range(0,size(lines)):
					at=np.append(at,np.double(timegm( strptime( lines[j][1:12]+'00', '%Y%m%d %H%M') )))
					adate=np.append(adate,date2num(datetime.datetime(time.gmtime(at[j])[0],time.gmtime(at[j])[1],time.gmtime(at[j])[2],time.gmtime(at[j])[3],time.gmtime(at[j])[4])))

					if len(lines[j])>0:
						ahs=np.append(ahs,np.float(lines[j][13:18]))
						ahspr=np.append(ahspr,np.float(lines[j][19:25]))
						atp=np.append(atp,np.float(lines[j][27:32]))
					else:
						ahs=np.append(ahs,np.nan)
						ahspr=np.append(ahspr,np.nan)
						atp=np.append(atp,np.nan)

				# build dictionary			
				result={'time':np.array(at).astype('double'),'date':np.array(adate).astype('double'),
				'station_name':np.array(stname),'hs':np.array(ahs),'hs_spr':np.array(ahspr),'tp':np.array(atp)}

				print("   Model data read, "+fname+", ts format.")
				return result
				del result,at,adate,ahs,ahspr,atp,tfile,lines

			elif 'glwu' in np.str(fname).split('/')[-1]:
				# great lakes ww3 format
				at=[];adate=[];ahs=[];al=[];atr=[];adir=[];aspr=[];atp=[];ap_dir=[];ap_spr=[]
				for j in range(0,size(lines)):
					at=np.append(at,np.double(timegm( strptime( lines[j][2:13]+'00', '%Y%m%d %H%M') )))
					adate=np.append(adate,date2num(datetime.datetime(time.gmtime(at[j])[0],time.gmtime(at[j])[1],time.gmtime(at[j])[2],time.gmtime(at[j])[3],time.gmtime(at[j])[4])))

					if len(lines[j])>0:
						ahs=np.append(ahs,np.float(lines[j][22:28]))
						al=np.append(al,np.float(lines[j][31:35]))
						atr=np.append(atr,np.float(lines[j][37:42]))
						adir=np.append(adir,np.float(lines[j][44:49]))
						aspr=np.append(aspr,np.float(lines[j][50:56]))
						atp=np.append(atp,np.float(lines[j][57:64]))
						ap_dir=np.append(ap_dir,np.float(lines[j][66:71]))
						ap_spr=np.append(ap_spr,np.float(lines[j][72:78]))
					else:
						ahs=np.append(ahs,np.nan)
						al=np.append(al,np.nan)
						atr=np.append(atr,np.nan)
						adir=np.append(adir,np.nan)
						aspr=np.append(aspr,np.nan)
						atp=np.append(atp,np.nan)
						ap_dir=np.append(ap_dir,np.nan)
						ap_spr=np.append(ap_spr,np.nan)

				atp[atp<0.01]=np.nan; atp=1./atp

				# build dictionary			
				result={'time':np.array(at).astype('double'),'date':np.array(adate).astype('double'),
				'station_name':np.array(stname),
				'hs':np.array(ahs),'l':np.array(al),
				'tm':np.array(atr),'dm':np.array(adir),
				'spr':np.array(aspr),'tp':np.array(atp),
				'dp':np.array(ap_dir),'peak_spr':np.array(ap_spr)}

				print("   Model data read, "+fname+", ts format.")
				return result
				del result,at,adate,ahs,al,atr,adir,aspr,atp,ap_dir,ap_spr,tfile,lines

	else:
		sys.exit(" Skipped file "+fname+" Not ts format.")


def station_tar(*args):
	'''
	WAVEWATCH III, station_tar operational point output, see https://www.ftp.ncep.noaa.gov/data/nccf/com/
	Input:  file name (example: gefs.wave.t00z.station_tar)
	Output: dictionary containing:
	  time(seconds since 1970),time(datetime64),station name; Arrays: hs, hs_spr, tp (gefs only)
	'''
	if len(args) == 1:
		fname=np.str(args[0])
	elif len(args) > 1:
		sys.exit(' Too many inputs')

	# confirm format
	if np.str(fname).split('/')[-1].split('.')[-1]=='station_tar':
		print("  reading ww3 station_tar file ...")
		import tarfile
		stname=[]

		try:
			tar = tarfile.open(fname)
		except:
			sys.exit('   Cannot open '+fname)
		else:
			for t in range(0,size(tar.getmembers())):
				# station names
				stname=np.append(stname,np.str(np.str(tar.getmembers()[t].name).split('/')[-1]).split('/')[-1].split('.')[-2])

				try:
					tfile=tar.extractfile(tar.getmembers()[t]); lines = tfile.readlines()[3::]
				except:
					print("   Cannot open "+tar.getmembers()[t].name)
				else:
					if t==0:
						# time array ----
						at=[]; adate=[]
						for j in range(0,size(lines)):
							at=np.append(at,np.double(timegm( strptime( np.str(lines[j])[3:14]+'00', '%Y%m%d %H%M') )))
							adate=np.append(adate,date2num(datetime.datetime(time.gmtime(at[j])[0],time.gmtime(at[j])[1],time.gmtime(at[j])[2],time.gmtime(at[j])[3],time.gmtime(at[j])[4])))

						# --------
						ahs=np.zeros((size(tar.getmembers()),at.shape[0]),'f')*np.nan
						ahspr=np.zeros((size(tar.getmembers()),at.shape[0]),'f')*np.nan
						atp=np.zeros((size(tar.getmembers()),at.shape[0]),'f')*np.nan

					auxhs=[]; auxhspr=[]; auxtp=[]
					for j in range(0,size(lines)):
						auxlines = np.str(lines[j]).replace("b'","")
						if len(lines[j])>0:
								auxhs=np.append(auxhs,np.float(auxlines[13:18]))
								auxhspr=np.append(auxhspr,np.float(auxlines[19:25]))
								auxtp=np.append(auxtp,np.float(auxlines[27:32]))
						else:
							auxhs=np.append(auxhs,np.nan)
							auxhspr=np.append(auxhspr,np.nan)
							auxtp=np.append(auxtp,np.nan)

					if ahs.shape[1]==auxhs.shape[0]:
						ahs[t,:]=np.array(auxhs)
						ahspr[t,:]=np.array(auxhspr)
						atp[t,:]=np.array(auxtp)
					else:
						print("   Time duration of "+tar.getmembers()[t]+" (in "+fname+") do not match the other stations. Mantained NaN.")

					del auxhs,auxhspr,auxtp,tfile,lines

			# build dictionary			
			result={'time':np.array(at).astype('double'),'date':np.array(adate).astype('double'),
			'station_name':np.array(stname),'hs':np.array(ahs),'hs_spr':np.array(ahspr),'tp':np.array(atp)}

			return result
			del result,tar,ahs,ahspr,atp,at,adate

			print(" Model data read, "+fname+", station_tar format.")

	else:
		sys.exit(" Skipped file "+fname+" Not station_tar format.")


# SPECTRA 

# Observations NDBC, netcdf format
def spec_ndbc(*args):
	'''
	Observations NDBC, wave spectrum, netcdf format
	Input: file name (example: 46047w2016.nc)
	Output: dictionary containing:
	  time(seconds since 1970),time(datetime64),lat,lon; Arrays: freq,dfreq,pspec,dmspec,dpspec,dirspec
	'''
	sk=1; deltatheta=np.int(10)
	if len(args) >= 1:
		fname=np.str(args[0])
	if len(args) >= 2:
		sk=np.int(args[1])
	if len(args) >= 3:
		deltatheta=np.int(args[3])
	if len(args) > 3:
		sys.exit(' Too many inputs')

	try:
		ds = xr.open_dataset(fname); f=nc.Dataset(fname)
	except:
		sys.exit(" Cannot open "+fname)
	else:
		btime = np.array(f.variables['time'][::sk]).astype('double')
		f.close(); del f
		bdate = ds['time'].values[::sk]
		blat = ds['latitude'].values[:]
		blon = ds['longitude'].values[:]
		freq = ds['frequency'].values[:]
		pspec = ds['spectral_wave_density'].values[::sk,:,0,0]
		dmspec = ds['mean_wave_dir'][::sk,:,0,0]
		dpspec = ds['principal_wave_dir'][::sk,:,0,0]	
		r1spec = ds['wave_spectrum_r1'][::sk,:,0,0]
		r2spec = ds['wave_spectrum_r2'][::sk,:,0,0]
		ds.close(); del ds
		# DF in frequency (dfreq), https://www.ndbc.noaa.gov/wavespectra.shtml
		dfreq=np.zeros(47,'f')
		dfreq[0]=0.010; dfreq[1:14]=0.005; dfreq[14:40]=0.010; dfreq[40::]=0.020
		pspec=np.array(pspec*dfreq)
		# Directional 2D Spectrum, https://www.ndbc.noaa.gov/measdes.shtml#swden , https://www.ndbc.noaa.gov/wavemeas.pdf
		theta = np.array(np.arange(0,360+0.1,deltatheta))
		# final directional wave spectrum (frequency X direction)
		dirspec = np.zeros((btime.shape[0],freq.shape[0],theta.shape[0]),'f')
		for t in range(0,btime.shape[0]):
			dirspec[t,:,:] = np.array([pspec[t,:]]).T * (1/pi)*(0.5+  np.array([r1spec[t,:]]).T * cos(np.array( np.array([theta])-np.array([dmspec[t,:]]).T )*(pi/180)) 
				+ np.array([r2spec[t,:]]).T*cos(2*np.array( np.array([theta]) - np.array([dpspec[t,:]]).T )*(pi/180)))

	# build dictionary			
	result={'time':btime,'date':bdate,'latitude':blat,'longitude':blon,
	'freq':freq,'deltafreq':dfreq,'pspec':pspec,'dmspec':dmspec,'dpspec':dpspec,
	'theta':theta,'dirspec':dirspec}

	return result
	del btime,bdate,blat,blon,freq,dfreq,pspec,dmspec,dpspec,theta,dirspec


# WAVEWATCH III spectra output, netcdf format
def spec_ww3(*args):
	'''
	WAVEWATCH III, wave spectrum, netcdf format
	Input: file name (example: ww3gefs.20160928_spec.nc), and station name (example: 41002)
	Output: dictionary containing:
	  time(seconds since 1970),time(datetime64),lat,lon; Arrays: freq,dfreq,pwst,d1sp,dire,dspec,wnds,wndd
	'''
	sk=1
	if len(args) < 2 :
		sys.exit(' Two inputs are required: file name and station name')
	if len(args) >= 2 :
		fname=np.str(args[0]); stname=np.str(args[1])
	if len(args) > 2 :
		sk=np.int(args[2])
	if len(args) > 3 :
		sys.exit(' Too many inputs')

	try:
		ds = xr.open_dataset(fname); f=nc.Dataset(fname)
	except:
		sys.exit(" Cannot open "+fname)
	else:

		mtime = np.array(f.variables['time'][::sk]*24*3600 + timegm( strptime(np.str(f.variables['time'].units).split(' ')[2][0:4]+'01010000', '%Y%m%d%H%M') )).astype('double')
		f.close(); del f

		auxstationname=ds['station_name'].values[:,:]; stationname=[]
		for i in range(0,auxstationname.shape[0]):
			stationname=np.append(stationname,"".join(np.array(auxstationname[i,:]).astype('str')))

		inds=np.where(stationname[:]==stname)
		if size(inds)>0:
			inds=np.int(inds[0][0]); stname=np.str(stationname[inds])
		else:
			sys.exit(' Station '+stname+' not included in the output file, or wrong station ID')

		# Spectrum
		dspec=np.array(ds['efth'][::sk,inds,:,:])
		# number of directions
		nd=dspec.shape[2]
		# number of frequencies
		nf=dspec.shape[1]
		# directions
		dire=np.array(ds['direction'].values[:])
		# frequencies
		freq=np.array(ds['frequency'].values[:])
		freq1=np.array(ds['frequency1'].values[:])
		freq2=np.array(ds['frequency2'].values[:])
		# DF in frequency (dfreq)
		dfreq=np.array(freq2 - freq1)
		# wind intensity and wind direction
		wnds=np.array(ds['wnd'].values[::sk,inds])
		wndd=np.array(ds['wnddir'].values[::sk,inds])
		# Time datetime64 array
		mdate=np.array(ds['time'].values[::sk])
		# water depth (constant in time)
		depth=np.nanmean(ds['dpt'].values[::sk,inds],axis=0)
		lon=np.array(np.nanmean(ds['longitude'].values[::sk,inds],axis=0))
		lat=np.array(np.nanmean(ds['latitude'].values[::sk,inds],axis=0))	
		
		ds.close(); del ds, auxstationname, inds, stationname
		# ------------------
		# 1D power spectrum
		pwst=np.zeros((dspec.shape[0],nf),'f')
		for t in range(0,dspec.shape[0]):
			for il in range(0,nf):
				pwst[t,il]=sum(dspec[t,il,:]*(2*np.pi)/nd)

			pwst[t,:]=pwst[t,:]*dfreq[:]

		# organizing directions  -----
		adspec=np.copy(dspec); inddire=int(np.where(dire==min(dire))[0][0])
		for t in range(0,dspec.shape[0]):
			adspec[t,:,0:nd-(inddire+1)]=dspec[t,:,(inddire+1):nd]
			adspec[t,:,nd-(inddire+1):nd]=dspec[t,:,0:(inddire+1)]
			for i in range(0,nd):
				dspec[t,:,i]=adspec[t,:,nd-i-1]

			adspec[t,:,0:int(nd/2)]=dspec[t,:,int(nd/2):nd]
			adspec[t,:,int(nd/2):nd]=dspec[t,:,0:int(nd/2)]
			dspec[t,:,:]=adspec[t,:,:]

		dire=np.sort(dire)

		# 1D directional spectrum
		d1sp=np.zeros((dspec.shape[0],nf),'f')
		for t in range(0,dspec.shape[0]):
			for il in range(0,nf):	
				a = np.sum(dspec[t,il,:] * np.array(np.sin((pi*dire)/180.)/np.sum(dspec[t,il,:])) )
				b = np.sum(dspec[t,il,:] * np.array(np.cos((pi*dire)/180.)/np.sum(dspec[t,il,:])) )
				aux = math.atan2(a,b)*(180./pi)
				if aux<0:
					aux=aux+360.
				
				d1sp[t,il]=np.float(aux)
				del a,b,aux

	# build dictionary			
	result={'time':mtime,'date':mdate,'latitude':lat,'longitude':lon,
	'wind_spd':wnds,'wind_dir':wndd,'freq':freq,'freq1':freq1,'freq2':freq2,
	'deltafreq':dfreq,'pspec':pwst,'theta':dire,'dmspec':d1sp,'dirspec':dspec}

	return result
	del mtime,mdate,lat,lon,wnds,wndd,freq,freq1,freq2,dfreq,pwst,dire,d1sp,dspec

