#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
wpanel2.py

VERSION AND LAST UPDATE:
 v1.0  06/03/2022

PURPOSE:
 Build a panel with multiple plots to compare two WAVEWATCHIII 
 similations, wrun1 and wrun2.
 In this panel it is shown:
   - Wind speed and wind current fields (common to both ww3 simulations);
   Hs (or any other field) for simulation 1 and 2, for comparison, and 
   the difference between sim1 and sim2.
   - A point comparison+evaluation is used, where a NDBC buoy and a ww3 
   point output is evaluated. Their 2D directional spectra and power
   spectra are presented.
   - Finaly a time-series plots of ww3 model and NDBC observations are
   included, for: Wind Speed, Hs, and Tp.

USAGE:
 variable of interest (wvar) must be informed below, default is Hs.
 paths of wave simulations have to be included
  and paths of buoy observations
 lat/lon limits/zoom-in for the plots can be edited

OUTPUT:
 png figures with multiple subplots: wave fiels, directional spectrum,
  power spectrum, and time-series of Wind Speed, Hs, and Tp.

DEPENDENCIES:
 See dependencies.py and the imports below.
 It uses function wread2.py

AUTHOR and DATE:
 06/03/2022: Ricardo M. Campos, first version.

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import matplotlib
matplotlib.use('Agg') # for backend plots, not for rendering in a window
import xarray as xr
import numpy as np
from pylab import *
from matplotlib.mlab import *
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import sys
import pandas as pd
import cartopy.crs as ccrs
import cartopy
from matplotlib import ticker
import matplotlib.colors as colors
from matplotlib.colors import BoundaryNorm
import warnings; warnings.filterwarnings("ignore")
import wread
palette = plt.cm.jet
dpalette = plt.cm.RdBu_r

# stname=np.str(sys.argv[1]) # stname="41047" 
# anh=np.str(sys.argv[2]) # anh=4.1 
stname="41049"
anh=4.1 # you can find anemometer height at the NDBC buoy's webpage
# lowest period (upper limit frequency) for the directional wave spectra (2D) polat plot
lper=4.5
# lat/lon limits/zoom-in for the plots
slat=np.array([-10,80]); slon=np.array([-168,-22])
# variable for the field plots
wvar="hs"
#   Paths of two ww3 simulations
# First WW3 simulation
wfile1="/data/ww3/gfs-d36.PR3/ww3.gfs-d36.PR3.glo_10mxt.20210924_20211024.nc"
wtab1="/data/ww3/gfs-d36.PR3/ww3.gfs-d36.tab.s2st2.20210924_20211024.nc"
wspec1="/data/ww3/gfs-d36.PR3/ww3.gfs-d36.spec.s1st3.20210924_20211024.nc"
trun1='d36.PR3'
# Second WW3 simulation
wfile2="/data/ww3/gfs-d48.PR3/ww3.gfs-d48.PR3.glo_10mxt.20210924_20211024.nc"
wtab2="/data/ww3/gfs-d48.PR3/ww3.gfs-d48.tab.s2st2.20210924_20211024.nc"
wspec2="/data/ww3/gfs-d48.PR3/ww3.gfs-d48.spec.s1st3.20210924_20211024.nc"
trun2='d48.PR3'
# NDBC data, time-series and spectrum
btab="/data/buoys/"+stname+"h2021.nc"
bspec="/data/buoys/"+stname+"w2021.nc"
# Skip time steps, to alleviate memory and reduce the number of consecutive plots.
sk=6

# READ DATA ********************

#   WW3 Fields ----------
# first simulation (main)
ds = xr.open_dataset(wfile1)
whs1 = ds[wvar]
wucur = ds['ucur']; wvcur = ds['vcur']
wuwnd = ds['uwnd']; wvwnd = ds['vwnd']
wtime = np.array(ds.time.values[:])
units_whs = np.str(ds[wvar].units)
units_wcur = np.str(ds['ucur'].units)
units_wwnd = np.str(ds['uwnd'].units)
lat = np.array(ds.latitude.values[:]); lon = np.array(ds.longitude.values[:])
ds.close(); del ds
# second simulation (main)
ds = xr.open_dataset(wfile2)
whs2 = ds[wvar]
ds.close(); del ds

#   TIME-SERIES
# Observations NDBC ----------------
tso_time,ignore,tso_lat,tso_lon,ignore,ignore,ignore,ignore,ignore,tso_wnd,ignore,tso_hs,tso_tm,tso_tp,tso_dm  = wread.tseriesnc_ndbc(btab)
# convert anemometer height to 10-meter level, approximation (DNVGL C-205).
tso_wnd = np.array( ((10./anh)**(0.1)) * tso_wnd)

# WAVEWATCH III ----------------
tsm_time,ignore,ignore,ignore,tsm_hs1,ignore,tsm_tr1,ignore,ignore,tsm_tp1,ignore,ignore = wread.tseriesnc_ww3b(wtab1,stname)
ignore,ignore,ignore,ignore,tsm_hs2,ignore,tsm_tr2,ignore,ignore,tsm_tp2,ignore,ignore = wread.tseriesnc_ww3b(wtab2,stname)

#   SPECTRA
# Observations NDBC ----------------------
spo_time,ignore,spo_lat,spo_lon,spo_freq,spo_dfreq,spo_pspec,spo_dmspec,ignore,spo_dire,spo_dspec = wread.spec_ndbc(bspec,1)
# WAVEWATCH III --------------------------
spm_time1,ignore,ignore,ignore,spm_freq,spm_freq1,spm_freq2,spm_dfreq,spm_pspec1,ignore,spm_dire1,spm_dspec1,tsm_wnd,ignore = wread.spec_ww3(wspec1,stname,sk)
spm_time2,ignore,ignore,ignore,ignore,ignore,ignore,ignore,spm_pspec2,ignore,spm_dire2,spm_dspec2,ignore,ignore = wread.spec_ww3(wspec2,stname,sk)

if np.array_equal(spm_time1,spm_time2)==True & np.array_equal(wtime,tsm_time)==True: 
	print(" Input data ok")
	spm_time=np.array(spm_time1); del spm_time1,spm_time2
else:
	sys.exit(' Error: two different model time-steps.')

# Time/index Intersect
ignore,ind1,ind2 = np.intersect1d(tsm_time,spm_time,assume_unique=False,return_indices=True)
#
ftime=np.array(tsm_time[ind1]); tsm_hs1=np.array(tsm_hs1[ind1]); tsm_hs2=np.array(tsm_hs2[ind1])
tsm_tr1=np.array(tsm_tr1[ind1]); tsm_tr2=np.array(tsm_tr2[ind1])
tsm_tp1=np.array(tsm_tp1[ind1]); tsm_tp2=np.array(tsm_tp2[ind1])
#
tsm_wnd=np.array(tsm_wnd[ind2])
spm_pspec1=np.array(spm_pspec1[ind2,:]); spm_pspec2=np.array(spm_pspec2[ind2,:])
spm_dspec1=np.array(spm_dspec1[ind2,:,:]); spm_dspec2=np.array(spm_dspec2[ind2,:,:])
# ********************

# FIGURES ----------------

# for the 2D polar plot:
indf=int(np.where(abs(spm_freq-(1/lper))==min(abs(spm_freq-(1/lper))))[0][0])
# WW3 model1
ndire1=np.zeros((spm_dire1.shape[0]+2),'f'); ndire1[1:-1]=spm_dire1[:]; ndire1[0]=0; ndire1[-1]=360
angle1 = np.radians(ndire1)
r1, theta1 = np.meshgrid(spm_freq[0:indf], angle1)
# WW3 model1
ndire2=np.zeros((spm_dire2.shape[0]+2),'f'); ndire2[1:-1]=spm_dire2[:]; ndire2[0]=0; ndire2[-1]=360
angle2 = np.radians(ndire2)
r2, theta2 = np.meshgrid(spm_freq[0:indf], angle2)
# NDBC buoy
ndire3=np.zeros((spo_dire.shape[0]+2),'f'); ndire3[1:-1]=spo_dire[:]; ndire3[0]=0; ndire3[-1]=360
angle3 = np.radians(ndire3)
indfb=int(np.where(abs(spo_freq-(1/lper))==min(abs(spo_freq-(1/lper))))[0][0])
r3, theta3 = np.meshgrid(spo_freq[0:indfb], angle3)
spo_freq1=np.array(spo_freq-spo_dfreq/2); spo_freq2=np.array(spo_freq+spo_dfreq/2)
spo_time = np.array(spo_time, dtype='datetime64[h]') # round to hour
# ----------------------------------------

levelshs = np.linspace(0.,9.,101); dlevelshs = np.linspace(-0.12,0.12,101)
levelswnd = np.linspace(0.,20.,101); levelscurr = np.linspace(0.,1.5,101)
print(" Ok. Starting plots ...")
for t in range(0,ftime.shape[0]):

	fig, axs = plt.subplots(nrows=3,ncols=4,subplot_kw={'projection': ccrs.PlateCarree()},figsize=(19,11))
	#   Wave Fields -----------------
	axs[0,1].set_extent([slon[0],slon[1],slat.min(),slat.max()], crs=ccrs.PlateCarree())  
	gl = axs[0,1].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
	gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
	gl.ylabels_right = False; gl.xlabels_top = False
	norm = BoundaryNorm(levelshs, ncolors=palette.N, clip=False)
	im = axs[0,1].pcolormesh(lon, lat, whs1[ind1[t],:,:],shading='flat',cmap=palette,norm=norm, zorder=2)
	axs[0,1].text(tso_lon, tso_lat,'*',color='white', size=12, zorder=3)
	axs[0,1].add_feature(cartopy.feature.OCEAN,facecolor=("white"))
	axs[0,1].add_feature(cartopy.feature.LAND,facecolor=("lightgrey"), edgecolor='grey',linewidth=0.5)
	axs[0,1].add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1)
	axs[0,1].coastlines(resolution='110m', color='grey',linewidth=0.5, linestyle='-', alpha=1)
	axs[0,1].set_title(wvar+' ('+units_whs+')    '+trun1+'    '+pd.to_datetime(ftime[t]).strftime('%Y/%m/%d %H')+'Z') 
	#
	axs[0,2].set_extent([slon[0],slon[1],slat.min(),slat.max()], crs=ccrs.PlateCarree())  
	gl = axs[0,2].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
	gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
	gl.ylabels_right = False; gl.xlabels_top = False
	im = axs[0,2].pcolormesh(lon, lat, whs2[ind1[t],:,:],shading='flat',cmap=palette,norm=norm, zorder=2)
	del norm
	axs[0,2].text(tso_lon, tso_lat,'*',color='white', size=12, zorder=3)
	axs[0,2].add_feature(cartopy.feature.OCEAN,facecolor=("white"))
	axs[0,2].add_feature(cartopy.feature.LAND,facecolor=("lightgrey"), edgecolor='grey',linewidth=0.5)
	axs[0,2].add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1)
	axs[0,2].coastlines(resolution='110m', color='grey',linewidth=0.5, linestyle='-', alpha=1)
	axs[0,2].set_title(wvar+' ('+units_whs+')    '+trun2+'    '+pd.to_datetime(ftime[t]).strftime('%Y/%m/%d %H')+'Z') 
	cax = axs[0,2].inset_axes([1.04, 0.2, 0.05, 0.6], transform=axs[0,2].transAxes)
	cbar = plt.colorbar(im, ax=axs[0,2], cax=cax, extend='max')
	tick_locator = ticker.MaxNLocator(nbins=6); cbar.locator = tick_locator; cbar.update_ticks()
	# Differences
	axs[0,3].set_extent([slon[0],slon[1],slat.min(),slat.max()], crs=ccrs.PlateCarree())  
	gl = axs[0,3].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
	gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
	gl.ylabels_right = False; gl.xlabels_top = False
	norm = BoundaryNorm(dlevelshs, ncolors=dpalette.N, clip=False)
	dim = axs[0,3].pcolormesh(lon, lat, np.array(whs1[ind1[t],:,:] - whs2[ind1[t],:,:]),shading='flat',cmap=dpalette,norm=norm, zorder=2)
	del norm
	axs[0,3].text(tso_lon, tso_lat,'*',color='k', size=12, zorder=3)
	axs[0,3].add_feature(cartopy.feature.OCEAN,facecolor=("white"))
	axs[0,3].add_feature(cartopy.feature.LAND,facecolor=("lightgrey"), edgecolor='grey',linewidth=0.5)
	axs[0,3].add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1)
	axs[0,3].coastlines(resolution='110m', color='grey',linewidth=0.5, linestyle='-', alpha=1)
	axs[0,3].set_title(wvar+' ('+units_whs+')  '+trun1+' - '+trun2) 
	cax = axs[0,3].inset_axes([1.04, 0.2, 0.05, 0.6], transform=axs[0,3].transAxes)
	cbar = plt.colorbar(dim, ax=axs[0,3], cax=cax, extend='both')
	tick_locator = ticker.MaxNLocator(nbins=6); cbar.locator = tick_locator; cbar.update_ticks()
	# Winds
	axs[0,0].set_extent([slon[0],slon[1],slat.min(),slat.max()], crs=ccrs.PlateCarree())  
	gl = axs[0,0].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
	gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
	gl.ylabels_right = False; gl.xlabels_top = False
	norm = BoundaryNorm(levelswnd, ncolors=palette.N, clip=False)
	im = axs[0,0].pcolormesh(lon, lat, np.sqrt(wuwnd[ind1[t],:,:]**2 + wvwnd[ind1[t],:,:]**2),shading='flat',cmap=palette,norm=norm, zorder=2)
	del norm
	axs[0,0].text(tso_lon, tso_lat,'*',color='white', size=12, zorder=3)
	axs[0,0].add_feature(cartopy.feature.OCEAN,facecolor=("white"))
	axs[0,0].add_feature(cartopy.feature.LAND,facecolor=("lightgrey"), edgecolor='grey',linewidth=0.5)
	axs[0,0].add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1)
	axs[0,0].coastlines(resolution='110m', color='grey',linewidth=0.5, linestyle='-', alpha=1)
	axs[0,0].set_title('10m Winds ('+units_wcur+')    '+pd.to_datetime(ftime[t]).strftime('%Y/%m/%d %H')+'Z') 
	cax = axs[0,0].inset_axes([1.04, 0.2, 0.05, 0.6], transform=axs[0,0].transAxes)
	cbar = plt.colorbar(im, ax=axs[0,0], cax=cax, extend='max')
	tick_locator = ticker.MaxNLocator(nbins=6); cbar.locator = tick_locator; cbar.update_ticks()
	# Currents
	axs[1,0].set_extent([slon[0],slon[1],slat.min(),slat.max()], crs=ccrs.PlateCarree())  
	gl = axs[1,0].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
	gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
	gl.ylabels_right = False; gl.xlabels_top = False
	norm = BoundaryNorm(levelscurr, ncolors=palette.N, clip=False)
	im = axs[1,0].pcolormesh(lon, lat, np.sqrt(wucur[ind1[t],:,:]**2 + wvcur[ind1[t],:,:]**2),shading='flat',cmap=palette,norm=norm, zorder=2)
	del norm
	axs[1,0].text(tso_lon, tso_lat,'*',color='k', size=12, zorder=3)
	axs[1,0].add_feature(cartopy.feature.OCEAN,facecolor=("white"))
	axs[1,0].add_feature(cartopy.feature.LAND,facecolor=("lightgrey"), edgecolor='grey',linewidth=0.5)
	axs[1,0].add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1)
	axs[1,0].coastlines(resolution='110m', color='grey',linewidth=0.5, linestyle='-', alpha=1)
	axs[1,0].set_title('Surface Currents ('+units_wcur+')    '+pd.to_datetime(ftime[t]).strftime('%Y/%m/%d %H')+'Z')
	cax = axs[1,0].inset_axes([1.04, 0.2, 0.05, 0.6], transform=axs[1,0].transAxes)
	cbar = plt.colorbar(im, ax=axs[1,0], cax=cax, extend='max')
	tick_locator = ticker.MaxNLocator(nbins=6); cbar.locator = tick_locator; cbar.update_ticks()
	# -----------------
	#   WW3 Wave Spectra -----------------
	# WW3 model1
	slevels = np.linspace(0.1,np.nanpercentile(np.c_[spm_dspec1[t,:,:],spm_dspec2[t,:,:]],99.99),201)
	ndspec1=np.zeros((spm_freq.shape[0],ndire1.shape[0]),'f')
	ndspec1[:,1:-1]=spm_dspec1[t,:,:]
	for i in range(0,spm_freq.shape[0]):
		ndspec1[i,-1]=float((ndspec1[i,-2]+ndspec1[i,1])/2.)
		ndspec1[i,0]=float((ndspec1[i,-2]+ndspec1[i,1])/2.)

	axs[1,1].remove()
	axs[1,1] = fig.add_subplot(3, 4, 6, projection='polar')
	axs[1,1].set_theta_zero_location('N')
	axs[1,1].set_theta_direction(-1)
	axs[1,1].set_rlabel_position(-135)
	axs[1,1].set_rticks([0.1,0.15,0.20]); axs[1,1].set_rmax(1/lper)
	im = axs[1,1].contourf(theta1, r1, ndspec1[0:indf,:].T,slevels,cmap=plt.cm.gist_stern_r,norm=colors.PowerNorm(gamma=0.5), extend="max")
	axs[1,1].set_title('WW3Spec '+stname+', '+trun1+' '+pd.to_datetime(ftime[t]).strftime('%Y/%m/%d %H')+'Z',size=11) 
	del im
	# WW3 model2
	ndspec2=np.zeros((spm_freq.shape[0],ndire2.shape[0]),'f')
	ndspec2[:,1:-1]=spm_dspec2[t,:,:]
	for i in range(0,spm_freq.shape[0]):
		ndspec2[i,-1]=float((ndspec2[i,-2]+ndspec2[i,1])/2.)
		ndspec2[i,0]=float((ndspec2[i,-2]+ndspec2[i,1])/2.)

	axs[1,2].remove()
	axs[1,2] = fig.add_subplot(3, 4, 7, projection='polar')
	axs[1,2].set_theta_zero_location('N')
	axs[1,2].set_theta_direction(-1)
	axs[1,2].set_rlabel_position(-135)
	axs[1,2].set_rticks([0.1,0.15,0.20]); axs[1,2].set_rmax(0.27)
	im = axs[1,2].contourf(theta2, r2, ndspec2[0:indf,:].T,slevels,cmap=plt.cm.gist_stern_r,norm=colors.PowerNorm(gamma=0.5), extend="max")
	axs[1,2].set_title('WW3Spec '+stname+', '+trun2+' '+pd.to_datetime(ftime[t]).strftime('%Y/%m/%d %H')+'Z',size=11) 
	#cax = axs[1,2].inset_axes([1.13, 0.2, 0.05, 0.6], transform=axs[1,2].transAxes)
	#cbar = plt.colorbar(im, ax=axs[1,2], cax=cax)
	#tick_locator = ticker.MaxNLocator(nbins=6); cbar.locator = tick_locator; cbar.update_ticks()
	del im
	# NDBC buoy
	axs[1,3].remove()
	indst=np.where( np.abs(spo_time-ftime[t]) == np.nanmin(np.abs(spo_time-ftime[t])) )
	if size(indst)>0:
		nslevels = np.linspace(0.001,np.nanpercentile(spo_dspec[indst[0],:,:],99.99),201)
		ndspec3=np.zeros((spo_freq.shape[0],ndire3.shape[0]),'f')
		ndspec3[:,1:-1]=spo_dspec[indst[0],:,:]
		for i in range(0,spo_freq.shape[0]):
			ndspec3[i,-1]=float((ndspec3[i,-2]+ndspec3[i,1])/2.)
			ndspec3[i,0]=float((ndspec3[i,-2]+ndspec3[i,1])/2.)

		axs[1,3] = fig.add_subplot(3, 4, 8, projection='polar')
		axs[1,3].set_theta_zero_location('N')
		axs[1,3].set_theta_direction(-1)
		axs[1,3].set_rlabel_position(-135)
		axs[1,3].set_rticks([0.1,0.15,0.20]); axs[1,3].set_rmax(0.27)
		im = axs[1,3].contourf(theta3, r3, ndspec3[0:indfb,:].T,nslevels,cmap=plt.cm.gist_stern_r,norm=colors.PowerNorm(gamma=0.5), extend="max")
		axs[1,3].set_title('NDBC '+stname+', '+pd.to_datetime(spo_time[indst[0]][0]).strftime('%Y/%m/%d %H')+'Z',size=11) 
		#cax = axs[1,3].inset_axes([1.13, 0.2, 0.05, 0.6], transform=axs[1,3].transAxes)
		#cbar = plt.colorbar(im, ax=axs[1,3], cax=cax)
		#tick_locator = ticker.MaxNLocator(nbins=6); cbar.locator = tick_locator; cbar.update_ticks()
		del im

		# *******   Spectral Interpolation Block  *******
		# Spectral Interpolation, Buoy and WW3, taking into account the spectral/frequency bands
		# Buoy
		spi_freq = np.arange(np.max([np.nanmin(spo_freq),np.nanmin(spm_freq)]), np.min([np.nanmax(spo_freq),np.nanmax(spm_freq)]), 0.0001) # very high resolution to guarantee the integral gives exactly the same origional Hs
		spoi_pspec = np.copy(spi_freq)*0.
		for i in range(1,spi_freq.shape[0]):		
			idf=np.where( ((spi_freq[i]-spo_freq1)>=0) & ((spi_freq[i]-spo_freq2)<=0) )[0][0]
			spoi_pspec[i] = 100.*np.nanmean(spo_pspec[indst[0],:],axis=0)[idf]*(spi_freq[i]-spi_freq[i-1])/spo_dfreq[idf]
			del idf

		# WW3
		spmi_pspec1 = np.copy(spi_freq)*0.
		spmi_pspec2 = np.copy(spi_freq)*0.
		for i in range(1,spi_freq.shape[0]):		
			idf=np.where( ((spi_freq[i]-spm_freq1)>=0) & ((spi_freq[i]-spm_freq2)<=0) )[0][0]
			spmi_pspec1[i] = 100.*spm_pspec1[t,idf]*(spi_freq[i]-spi_freq[i-1])/spm_dfreq[idf]
			spmi_pspec2[i] = 100.*spm_pspec2[t,idf]*(spi_freq[i]-spi_freq[i-1])/spm_dfreq[idf]
			del idf

		# *********************
		axs[2,3].remove()
		axs[2,3] = fig.add_subplot(3, 4, 12, projection='rectilinear')
		axs[2,3].plot(spi_freq[1:3650],spmi_pspec1[1:3650], color='b', linestyle='--',linewidth=2.0, label=trun1, zorder=3)
		axs[2,3].plot(spi_freq[1:3650],spmi_pspec2[1:3650], color='r', linestyle='-.',linewidth=2.0, label=trun2, zorder=3)
		axs[2,3].plot(spi_freq[1:3650],spoi_pspec[1:3650], color='grey', linestyle='-',linewidth=0.5, label='buoy', zorder=2)
		axs[2,3].legend(loc='best')
		axs[2,3].fill_between(spi_freq[1:3650], 0.,spoi_pspec[1:3650], color='silver', alpha=0.7,label='buoy', zorder=1)
		axs[2,3].grid(linewidth=0.5, color='grey', alpha=0.5, linestyle='--', zorder=1)
		axs[2,3].set_xlabel('Frequency (Hz)', fontsize=12); axs[2,3].set_ylabel('Power Spectrum (m$^2$/Hz)', fontsize=12) 
		axs[2,3].set_ylim(ymin = -0.001)
		axs[2,3].set_title('PowerSpec '+stname+' '+pd.to_datetime(ftime[t]).strftime('%Y/%m/%d %H')+'Z')
		axs[2,3].axis('tight')

	# -----------------
	# Time-Series -----------------
	axs[2,0].remove(); axs[2,1].remove(); axs[2,2].remove()
	iaux=np.intersect1d(ftime[np.max([0,t-12]):np.min([np.int(ftime.shape[0]-1),t+13])], np.array(tso_time, dtype='datetime64[h]'), assume_unique=False, return_indices=True)
	axs[2,0] = fig.add_subplot(3, 4, 9, projection='rectilinear')
	axs[2,0].plot_date(ftime[np.max([0,t-12]):np.min([np.int(ftime.shape[0]-1),t+13])],tsm_wnd[np.max([0,t-12]):np.min([np.int(ftime.shape[0]-1),t+13])],color='navy', linestyle='--',marker='',linewidth=2.0, label='GFS winds', zorder=3)
	if size(iaux)>0:
		axs[2,0].fill_between(tso_time[np.min(iaux[2]):np.max(iaux[2])+1], 0., tso_wnd[np.min(iaux[2]):np.max(iaux[2])+1], color='silver',alpha=0.5,zorder=1)
		axs[2,0].plot_date(tso_time[np.min(iaux[2]):np.max(iaux[2])+1],tso_wnd[np.min(iaux[2]):np.max(iaux[2])+1],color='dimgrey',marker='.',label='buoy', zorder=2)

	axs[2,0].xaxis.set_major_formatter( DateFormatter('%b%d') ); axs[2,0].fmt_xdata = DateFormatter('%b%d')
	axs[2,0].axvline(x=ftime[t],ls='--',color='grey',linewidth=1.,alpha=0.9,zorder=1)
	axs[2,0].grid(linewidth=0.5, color='grey', alpha=0.5, linestyle='--', zorder=1)
	axs[2,0].xaxis.set_major_locator(mdates.DayLocator())
	axs[2,0].legend(loc='best', fontsize=9)
	axs[2,0].set_xlabel('Date', fontsize=12); axs[2,0].set_ylabel('Wind Speed (m/s)', fontsize=12) 
	axs[2,0].axis('tight')
	axs[2,0].set_xlim( ftime[np.max([0,t-12]):np.min([np.int(ftime.shape[0]-1),t+13])].min(), ftime[np.max([0,t-12]):np.min([np.int(ftime.shape[0]-1),t+13])].max() )
	#
	axs[2,1] = fig.add_subplot(3, 4, 10, projection='rectilinear')
	axs[2,1].plot_date(ftime[np.max([0,t-12]):np.min([np.int(ftime.shape[0]-1),t+13])],tsm_hs1[np.max([0,t-12]):np.min([np.int(ftime.shape[0]-1),t+13])],color='b', linestyle='--',marker='',linewidth=2.0, label=trun1, zorder=3)
	axs[2,1].plot_date(ftime[np.max([0,t-12]):np.min([np.int(ftime.shape[0]-1),t+13])],tsm_hs2[np.max([0,t-12]):np.min([np.int(ftime.shape[0]-1),t+13])],color='r', linestyle='-.',marker='',linewidth=2.0, label=trun2, zorder=3)
	if size(iaux)>0:
		axs[2,1].fill_between(tso_time[np.min(iaux[2]):np.max(iaux[2])+1], 0., tso_hs[np.min(iaux[2]):np.max(iaux[2])+1], color='silver',alpha=0.5,zorder=1)
		axs[2,1].plot_date(tso_time[np.min(iaux[2]):np.max(iaux[2])+1],tso_hs[np.min(iaux[2]):np.max(iaux[2])+1],color='dimgrey',marker='.',label='buoy', zorder=2)
	axs[2,1].xaxis.set_major_formatter( DateFormatter('%b%d') ); axs[2,1].fmt_xdata = DateFormatter('%b%d')
	axs[2,1].axvline(x=ftime[t],ls='--',color='grey',linewidth=1.,alpha=0.9,zorder=1)
	axs[2,1].grid(linewidth=0.5, color='grey', alpha=0.5, linestyle='--', zorder=1)
	axs[2,1].xaxis.set_major_locator(mdates.DayLocator())
	axs[2,1].legend(loc='best', fontsize=9)
	axs[2,1].set_xlabel('Date', fontsize=12); axs[2,1].set_ylabel('Hs (m)', fontsize=12) 
	axs[2,1].axis('tight')
	axs[2,1].set_xlim( ftime[np.max([0,t-12]):np.min([np.int(ftime.shape[0]-1),t+13])].min(), ftime[np.max([0,t-12]):np.min([np.int(ftime.shape[0]-1),t+13])].max() )
	#
	axs[2,2] = fig.add_subplot(3, 4, 11, projection='rectilinear')
	axs[2,2].plot_date(ftime[np.max([0,t-12]):np.min([np.int(ftime.shape[0]-1),t+13])],tsm_tp1[np.max([0,t-12]):np.min([np.int(ftime.shape[0]-1),t+13])],color='b', linestyle='--',marker='',linewidth=2.0, label=trun1, zorder=3)
	axs[2,2].plot_date(ftime[np.max([0,t-12]):np.min([np.int(ftime.shape[0]-1),t+13])],tsm_tp2[np.max([0,t-12]):np.min([np.int(ftime.shape[0]-1),t+13])],color='r', linestyle='-.',marker='',linewidth=2.0, label=trun2, zorder=3)
	if size(iaux)>0:
		axs[2,2].fill_between(tso_time[np.min(iaux[2]):np.max(iaux[2])+1], 0., tso_tp[np.min(iaux[2]):np.max(iaux[2])+1], color='silver',alpha=0.5,zorder=1)
		axs[2,2].plot_date(tso_time[np.min(iaux[2]):np.max(iaux[2])+1],tso_tp[np.min(iaux[2]):np.max(iaux[2])+1],color='dimgrey',marker='.',label='buoy', zorder=2)
	axs[2,2].xaxis.set_major_formatter( DateFormatter('%b%d') ); axs[2,2].fmt_xdata = DateFormatter('%b%d')
	axs[2,2].axvline(x=ftime[t],ls='--',color='grey',linewidth=1.,alpha=0.9,zorder=1)
	axs[2,2].grid(linewidth=0.5, color='grey', alpha=0.5, linestyle='--', zorder=1)
	axs[2,2].xaxis.set_major_locator(mdates.DayLocator())
	axs[2,2].legend(loc='best', fontsize=9)
	axs[2,2].set_xlabel('Date', fontsize=12); axs[2,2].set_ylabel('Tp (s)', fontsize=12) 
	axs[2,2].axis('tight')
	axs[2,2].set_xlim( ftime[np.max([0,t-12]):np.min([np.int(ftime.shape[0]-1),t+13])].min(), ftime[np.max([0,t-12]):np.min([np.int(ftime.shape[0]-1),t+13])].max() )

	fig.canvas.draw() # https://github.com/SciTools/cartopy/issues/1207
	fig.tight_layout()
	plt.savefig('WPanel_'+stname+'_'+np.str(pd.to_datetime(ftime[t]).strftime('%Y%m%d%H'))+'_'+trun1+'_'+trun2+'.png', dpi=300, facecolor='w', edgecolor='w',
		orientation='portrait', papertype=None, format='png',transparent=False, pad_inches=0.1)

	plt.close('all'); del axs, fig, indst, iaux

print(" Done.")
