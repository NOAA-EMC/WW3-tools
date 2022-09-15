#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
wread.py

VERSION AND LAST UPDATE:
 v1.0  04/04/2022

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
   tseriesnc_ndbc
   tseriesnc_copernicus
   tseriesnc_ww3
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


# TIME-SERIES

# Observations NDBC, netcdf format
def tseriesnc_ndbc(*args):
	'''
	Observations NDBC, time series/table, netcdf format
	Input: file name (example: 46047h2016.nc)
	Output: Values: time(datetime64),time(seconds since 1970),lat,lon; Arrays: sst,mslp,dwp,tmp,gst,wsp,wdir,hs,tm,tp,dm
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
		bdate = ds['time'].values[:]
		blat = ds['latitude'].values[:]
		blon = ds['longitude'].values[:]
		bsst = ds['sea_surface_temperature'].values[:,0,0]
		bmslp = ds['air_pressure'].values[:,0,0]
		bdwp = ds['dewpt_temperature'].values[:,0,0]
		btmp = ds['air_temperature'].values[:,0,0]
		bgst = ds['gust'].values[:,0,0]	
		bwsp = ds['wind_spd'].values[:,0,0]
		bwdir = ds['wind_dir'].values[:]
		bhs = ds['wave_height'].values[:,0,0]
		# btm = ds['average_wpd'].values[:,0,0]
		# btp = ds['dominant_wpd'].values[:,0,0]
		bdm = ds['mean_wave_dir'].values[:,0,0]
		ds.close(); del ds

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

		return bdate,btime,blat,blon,bsst,bmslp,bdwp,btmp,bgst,bwsp,bwdir,bhs,btm,btp,bdm

# Observations Copernicus, netcdf format
def tseriesnc_copernicus(*args):
	'''
	Observations NDBC, time series/table, netcdf format
	Input: file name (example: 46047h2016.nc)
	Output: Values: time(datetime64),time(seconds since 1970),lat,lon; Arrays: sst,mslp,dwp,tmp,gst,wsp,wdir,hs,tm,tp,dm
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
		bdate = ds['TIME'].values[:]
		blat = np.nanmean(ds['LATITUDE'].values[:])
		blon = np.nanmean(ds['LONGITUDE'].values[:])
		bsst = np.nanmean(ds['TEMP'].values[:,:],axis=1) # SST 
		bmslp = np.nanmean(ds['ATMS'].values[:,:],axis=1) # Pressure
		bdwp = np.nanmean(ds['DEWT'].values[:,:],axis=1) # Dewpoint
		btmp = np.nanmean(ds['DRYT'].values[:,:],axis=1) # temperature
		bgst = np.nanmean(ds['GSPD'].values[:,:],axis=1) # gust
		bwsp = np.nanmean(ds['WSPD'].values[:,:],axis=1) # wind speed		
		bwdir = np.nanmean(ds['WDIR'].values[:,:],axis=1) # wind direction		
		bhs = np.nanmean(ds['VHM0'].values[:,:],axis=1) # Hs
		btm = np.nanmean(ds['VTM02'].values[:,:],axis=1) # Tm
		btp = np.nanmean(ds['VTPK'].values[:,:],axis=1) # Tp
		bdm = np.nanmean(ds['VMDR'].values[:,:],axis=1) # Mean direction

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
		ds.close(); del ds

		return bdate,btime,blat,blon,bsst,bmslp,bdwp,btmp,bgst,bwsp,bwdir,bhs,btm,btp,bdm

# WAVEWATCH III point output, netcdf format
def tseriesnc_ww3(*args):
	'''
	WAVEWATCH III, time series/table, netcdf format
	Input:  file name (example: ww3gefs.20160928_tab.nc), and station name (example: 41002)
		  point output file created with ww3_ounf last lines options:
		  T 1 
		  2
		  0
		  T
		  2
	Output: Values: time(datetime64),time(seconds since 1970),lat,lon; Arrays: hs, tm, tp, dm, dp, spr, lm
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

		mdate = ds['time'].values[:]
		mlat = np.nanmean(ds['latitude'].values[:,inds])
		mlon = np.nanmean(ds['longitude'].values[:,inds])
		mhs = ds['hs'].values[:,inds]
		mtp = np.zeros(mhs.shape[0],'f')*np.nan
		mtp[ds['fp'].values[:,inds]>0.0] = 1./ds['fp'].values[:,inds][ds['fp'].values[:,inds]>0.0]
		mtm = ds['tr'].values[:,inds]
		mdp = ds['th1p'].values[:,inds]
		mdm = ds['th1m'].values[:,inds]
		spr = ds['sth1m'].values[:,inds]
		lm = ds['lm'].values[:,inds]
		ds.close(); del ds

	return mdate,mtime,mlat,mlon,mhs,mtm,mtp,mdm,mdp,spr,lm

def tseriesnc_ww3b(*args):
	'''
	WAVEWATCH III, time series/table, netcdf format
	Input:  file name (example: ww3gefs.20160928_tab.nc), and station name (example: 41002)
		  point output file created with ww3_ounf last lines options:

	Output: Values: time(datetime64),time(seconds since 1970),lat,lon; Arrays: hs, lm, tr, th1p, sth1p, tp, th1m, sth1m
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

		mdate = ds['time'].values[:]
		mlat = np.nanmean(ds['latitude'].values[:,inds])
		mlon = np.nanmean(ds['longitude'].values[:,inds])
		mhs = ds['hs'].values[:,inds]
		mlm = ds['lm'].values[:,inds]
		mtr = ds['tr'].values[:,inds]
		mth1p = ds['th1p'].values[:,inds]
		msth1p = ds['sth1p'].values[:,inds]
		mtp = np.zeros(mhs.shape[0],'f')*np.nan
		mtp[ds['fp'].values[:,inds]>0.0] = 1./ds['fp'].values[:,inds][ds['fp'].values[:,inds]>0.0]
		mth1m = ds['th1m'].values[:,inds]
		msth1m = ds['sth1m'].values[:,inds]
		ds.close(); del ds

	return mdate,mtime,mlat,mlon,mhs,mlm,mtr,mth1p,msth1p,mtp,mth1m,msth1m

# SPECTRA 

# Observations NDBC, netcdf format
def spec_ndbc(*args):
	'''
	Observations NDBC, wave spectrum, netcdf format
	Input: file name (example: 46047w2016.nc)
	Output: Values: time(datetime64),time(seconds since 1970),lat,lon; Arrays: freq,dfreq,pspec,dmspec,dpspec,dirspec
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

	return bdate,btime,blat,blon,freq,dfreq,pspec,dmspec,dpspec,theta,dirspec


# WAVEWATCH III spectra output, netcdf format
def spec_ww3(*args):
	'''
	WAVEWATCH III, wave spectrum, netcdf format
	Input: file name (example: ww3gefs.20160928_spec.nc), and station name (example: 41002)
	Output: Values: time(datetime64),time(seconds since 1970),lat,lon; Arrays: freq,dfreq,pwst,d1sp,dire,dspec,wnds,wndd
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

	return mdate,mtime,lat,lon,freq,freq1,freq2,dfreq,pwst,d1sp,dire,dspec,wnds,wndd


