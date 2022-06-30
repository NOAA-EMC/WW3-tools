#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
wpanel1.py

VERSION AND LAST UPDATE:
 v1.0  04/04/2022

PURPOSE:
 Build a panel with multiple plots to compare two WAVEWATCHIII 
 similations, wrun1 and wrun2.
 Initial and final date time is defined by dt (see below) a 
 domain-selection(zoom-in) is done for the wave field plots by 
 defining slat and slon

USAGE:
 variable of interest (wvar) must be informed below
 paths of wave simulations have to be included (wrun)
  and paths of buoy observations
 lat/lon limits/zoom-in for the plots can be edited

OUTPUT:
 png figures with multiple subplots: wave fiels, directional spectrum,
  power spectrum, and time-series of Hs, Tm, and Dm.

DEPENDENCIES:
 See dependencies.py and the imports below.
 It uses function wread.py

AUTHOR and DATE:
 04/04/2022: Ricardo M. Campos, first version.

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import matplotlib
matplotlib.use('Agg')  # uncomment this for backend plots, not for rendering in a window
import xarray as xr
import numpy as np
from pylab import *
from matplotlib.mlab import *
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import sys
import pandas as pd
import cartopy.crs as ccrs
import cartopy
from scipy.ndimage.filters import gaussian_filter
from matplotlib import ticker
import matplotlib.colors as colors
from matplotlib.colors import BoundaryNorm
# import pickle
import warnings; warnings.filterwarnings("ignore")
import wread
palette = plt.cm.jet
dpalette = plt.cm.RdBu_r

stname=np.str(sys.argv[1]) # stname="41047" 

# lowest period (upper limit frequency) for the directional wave spectra (2D) polat plot
lper=4.5
# lat/lon limits/zoom-in for the plots
slat=np.array([5,65]); slon=np.array([-90,-20])
# variable for the field plots
wvar="hs"
# Paths
wrun1="/data/ww3/c00"; trun1="PR3.UQ.WDTH1.5" # Title of run
wrun2="/data/ww3/c00"; trun2="PR3.UQ.WDTH1.5_48D"
pbuoys="/data/buoys"

# READ DATA ********************
dt = np.arange(datetime(2016,10,4), datetime(2016,11,3), timedelta(days=1)).astype(datetime) # initial/final date of hindcast (each file with 1 day)
it=0
for t in range(0,size(dt)):

	cdate = pd.to_datetime(dt[t]).strftime('%Y%m%d')

	#   WW3 Fields ----------
	fname=wrun1+"/ww3gefs."+cdate+"_field.nc"
	ds = xr.open_dataset(fname)
	if it==0:
		wdata1 = np.array(ds[wvar].values[:,:,:])
		wtime = np.array(ds.time.values)
		units_wdata = np.str(ds[wvar].units)
		lat = np.array(ds.latitude.values[:]); lon = np.array(ds.longitude.values[:])
	else:
		wdata1 = np.append(wdata1,np.array(ds[wvar].values[:,:,:]),axis=0)
		wtime = np.append(wtime,np.array(ds.time.values))

	ds.close(); del ds
	# 
	fname=wrun2+"/ww3gefs."+cdate+"_field.nc"
	ds = xr.open_dataset(fname)
	if it==0:
		wdata2 = np.array(ds[wvar].values[:,:,:])
	else:
		wdata2 = np.append(wdata2,np.array(ds[wvar].values[:,:,:]),axis=0)

	ds.close(); del ds


	#   TIME-SERIES
	# Observations NDBC ----------------
	if it==0:
		fname=np.str(pbuoys+"/"+stname+"h2016.nc")
		tso_time,lixo,tso_lat,tso_lon,lixo,lixo,lixo,lixo,lixo,tso_wnd,lixo,tso_hs,tso_tm,lixo,tso_dm  = wread.tseriesnc_ndbc(fname)
		# round to hour, remove repeated records, running_average of 3h
		# tso_time = np.array(tso_time, dtype='datetime64[h]') # round to hour
		del fname


	# WAVEWATCH III ----------------
	if it==0:
		fname=np.str(wrun1+"/ww3gefs."+cdate+"_tab.nc")
		tsm_time,lixo,lixo,lixo,tsm_hs1,tsm_tm1,lixo,tsm_dm1,lixo,lixo,lixo = wread.tseriesnc_ww3(fname,stname)
		del fname
	else:
		fname=np.str(wrun1+"/ww3gefs."+cdate+"_tab.nc")
		atsm_time,lixo,lixo,lixo,atsm_hs1,atsm_tm1,lixo,atsm_dm1,lixo,lixo,lixo = wread.tseriesnc_ww3(fname,stname)
		del fname
		tsm_time = np.append(tsm_time,atsm_time)
		tsm_hs1 = np.append(tsm_hs1,atsm_hs1)
		tsm_tm1 = np.append(tsm_tm1,atsm_tm1)
		tsm_dm1 = np.append(tsm_dm1,atsm_dm1)

	if it==0:
		fname=np.str(wrun2+"/ww3gefs."+cdate+"_tab.nc")
		lixo,lixo,lixo,lixo,tsm_hs2,tsm_tm2,lixo,tsm_dm2,lixo,lixo,lixo = wread.tseriesnc_ww3(fname,stname)
		del fname
	else:
		fname=np.str(wrun2+"/ww3gefs."+cdate+"_tab.nc")
		lixo,lixo,lixo,lixo,atsm_hs2,atsm_tm2,lixo,atsm_dm2,lixo,lixo,lixo = wread.tseriesnc_ww3(fname,stname)
		del fname
		tsm_hs2 = np.append(tsm_hs2,atsm_hs2)
		tsm_tm2 = np.append(tsm_tm2,atsm_tm2)
		tsm_dm2 = np.append(tsm_dm2,atsm_dm2)

	# ---------------------------------

	# SPECTRA
	# Observations NDBC ----------------------
	if it==0:
		fname = pbuoys+"/"+stname+"w2016.nc"
		spo_time,lixo,spo_lat,spo_lon,spo_freq,spo_dfreq,spo_pspec,spo_dmspec,lixo,lixo,lixo = wread.spec_ndbc(fname)
		spo_freq1=np.array(spo_freq-spo_dfreq/2); spo_freq2=np.array(spo_freq+spo_dfreq/2)
		spo_time = np.array(spo_time, dtype='datetime64[h]') # round to hour
		del fname
		indfb=int(np.where(abs(spo_freq-(1/lper))==min(abs(spo_freq-(1/lper))))[0][0])

	# ----------------------------------------

	# WAVEWATCH III --------------------------
	if it==0:
		fname=wrun1+"/ww3gefs."+cdate+"_spec.nc"
		spm_time,lixo,lixo,lixo,spm_freq,spm_freq1,spm_freq2,spm_dfreq,spm_pspec1,spm_dmspec1,spm_dire,spm_dspec1,lixo,lixo = wread.spec_ww3(fname,stname)
		del fname
	else:
		fname=wrun1+"/ww3gefs."+cdate+"_spec.nc"
		aspm_time,lixo,lixo,lixo,lixo,lixo,lixo,lixo,aspm_pspec1,aspm_dmspec1,lixo,aspm_dspec1,lixo,lixo = wread.spec_ww3(fname,stname)
		del fname
		spm_time = np.append(spm_time,aspm_time)
		spm_pspec1 = np.append(spm_pspec1,aspm_pspec1,axis=0)
		spm_dmspec1 = np.append(spm_dmspec1,aspm_dmspec1,axis=0)
		spm_dspec1 = np.append(spm_dspec1,aspm_dspec1,axis=0)

	if it==0:
		fname=wrun2+"/ww3gefs."+cdate+"_spec.nc"
		lixo,lixo,lixo,lixo,lixo,lixo,lixo,lixo,spm_pspec2,spm_dmspec2,lixo,spm_dspec2,lixo,lixo = wread.spec_ww3(fname,stname)
		del fname
	else:
		fname=wrun2+"/ww3gefs."+cdate+"_spec.nc"
		lixo,lixo,lixo,lixo,lixo,lixo,lixo,lixo,aspm_pspec2,aspm_dmspec2,lixo,aspm_dspec2,lixo,lixo = wread.spec_ww3(fname,stname)
		del fname
		spm_pspec2 = np.append(spm_pspec2,aspm_pspec2,axis=0)
		spm_dmspec2 = np.append(spm_dmspec2,aspm_dmspec2,axis=0)
		spm_dspec2 = np.append(spm_dspec2,aspm_dspec2,axis=0)

	it=it+1

# FIGURES ----------------

# for the 2D polar plot:
indf=int(np.where(abs(spm_freq-(1/lper))==min(abs(spm_freq-(1/lper))))[0][0])
ndire=np.zeros((spm_dire.shape[0]+2),'f'); ndire[1:-1]=spm_dire[:]; ndire[0]=0; ndire[-1]=360
angle = np.radians(ndire)
r, theta = np.meshgrid(spm_freq[0:indf], angle)
# ----------------------------------------

levels = np.linspace(np.nanmin(np.c_[wdata1,wdata2]),np.nanpercentile(np.c_[wdata1,wdata2],99.9),101)
dlevels = np.linspace(-0.2,0.2,101); dlevelsp = np.linspace(0.,10.,101)

for t in range(20,size(wtime)-20):

	fig, axs = plt.subplots(nrows=3,ncols=4,subplot_kw={'projection': ccrs.PlateCarree()},figsize=(19,11))
	# axs=axs.flatten()
	# Wave Fields -----------------
	axs[0,0].set_extent([slon[0],slon[1],slat.min(),slat.max()], crs=ccrs.PlateCarree())  
	gl = axs[0,0].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
	gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
	gl.ylabels_right = False; gl.xlabels_top = False
	# im = axs[0,0].contourf(lon,lat,wdata1[t,:,:],levels,cmap=palette,extend="max", zorder=2)
	norm = BoundaryNorm(levels, ncolors=palette.N, clip=False)
	im = axs[0,0].pcolormesh(lon, lat, wdata1[t,:,:],shading='flat',cmap=palette,norm=norm, zorder=2)
	axs[0,0].text(tso_lon, tso_lat,'*',color='k', size=12, zorder=3)
	axs[0,0].add_feature(cartopy.feature.OCEAN,facecolor=("white"))
	axs[0,0].add_feature(cartopy.feature.LAND,facecolor=("lightgrey"), edgecolor='grey',linewidth=0.5)
	axs[0,0].add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1)
	axs[0,0].coastlines(resolution='110m', color='grey',linewidth=0.5, linestyle='-', alpha=1)
	axs[0,0].set_title(wvar+' ('+units_wdata+')    '+trun1+'    '+pd.to_datetime(wtime[:][t]).strftime('%Y/%m/%d %H')+'Z') 
	#
	axs[0,1].set_extent([slon[0],slon[1],slat.min(),slat.max()], crs=ccrs.PlateCarree())  
	gl = axs[0,1].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
	gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
	gl.ylabels_right = False; gl.xlabels_top = False
	# im = axs[0,1].contourf(lon,lat,wdata2[t,:,:],levels,cmap=palette,extend="max", zorder=2)
	im = axs[0,1].pcolormesh(lon, lat, wdata2[t,:,:],shading='flat',cmap=palette,norm=norm, zorder=2)
	del norm
	axs[0,1].text(tso_lon, tso_lat,'*',color='k', size=12, zorder=3)
	axs[0,1].add_feature(cartopy.feature.OCEAN,facecolor=("white"))
	axs[0,1].add_feature(cartopy.feature.LAND,facecolor=("lightgrey"), edgecolor='grey',linewidth=0.5)
	axs[0,1].add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1)
	axs[0,1].coastlines(resolution='110m', color='grey',linewidth=0.5, linestyle='-', alpha=1)
	axs[0,1].set_title(wvar+' ('+units_wdata+')    '+trun2+'    '+pd.to_datetime(wtime[:][t]).strftime('%Y/%m/%d %H')+'Z') 
	cax = axs[0,1].inset_axes([1.04, 0.2, 0.05, 0.6], transform=axs[0,1].transAxes)
	cbar = plt.colorbar(im, ax=axs[0,1], cax=cax, extend='max')
	tick_locator = ticker.MaxNLocator(nbins=6); cbar.locator = tick_locator; cbar.update_ticks()
	#
	axs[0,2].set_extent([slon[0],slon[1],slat.min(),slat.max()], crs=ccrs.PlateCarree())  
	gl = axs[0,2].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
	gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
	gl.ylabels_right = False; gl.xlabels_top = False
	# dim = axs[0,2].contourf(lon,lat,np.array(wdata1 - wdata2)[t,:,:],dlevels,cmap=dpalette,extend="both", zorder=2)
	norm = BoundaryNorm(dlevels, ncolors=dpalette.N, clip=False)
	dim = axs[0,2].pcolormesh(lon, lat, np.array(wdata1 - wdata2)[t,:,:],shading='flat',cmap=dpalette,norm=norm, zorder=2)
	del norm
	axs[0,2].text(tso_lon, tso_lat,'*',color='k', size=12, zorder=3)
	axs[0,2].add_feature(cartopy.feature.OCEAN,facecolor=("white"))
	axs[0,2].add_feature(cartopy.feature.LAND,facecolor=("lightgrey"), edgecolor='grey',linewidth=0.5)
	axs[0,2].add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1)
	axs[0,2].coastlines(resolution='110m', color='grey',linewidth=0.5, linestyle='-', alpha=1)
	axs[0,2].set_title(wvar+' ('+units_wdata+')  '+trun1+' - '+trun2) 
	cax = axs[0,2].inset_axes([1.04, 0.2, 0.05, 0.6], transform=axs[0,2].transAxes)
	cbar = plt.colorbar(dim, ax=axs[0,2], cax=cax, extend='both')
	tick_locator = ticker.MaxNLocator(nbins=6); cbar.locator = tick_locator; cbar.update_ticks()
	#
	axs[0,3].remove()
	fig.text(0.35,1.3,"hs_max "+trun1+": "+np.str(np.round(np.nanmax(wdata1[t,:,:]),2)),color='k', size=14, horizontalalignment='left', verticalalignment='center', transform=axs[0,3].transAxes)
	fig.text(0.35,1.1,"hs_max "+trun2+": "+np.str(np.round(np.nanmax(wdata2[t,:,:]),2)),color='k', size=14, horizontalalignment='left', verticalalignment='center', transform=axs[0,3].transAxes)
	fig.text(0.35,0.9,"hs_mean "+trun1+": "+np.str(np.round(np.nanmean(wdata1[t,:,:]),2)),color='k', size=14, horizontalalignment='left', verticalalignment='center', transform=axs[0,3].transAxes)
	fig.text(0.35,0.7,"hs_mean "+trun2+": "+np.str(np.round(np.nanmean(wdata2[t,:,:]),2)),color='k', size=14, horizontalalignment='left', verticalalignment='center', transform=axs[0,3].transAxes)
	# -----------------
	# WW3 Wave Spectra -----------------
	slevels = np.linspace(0.1,np.nanpercentile(np.c_[spm_dspec1[t,:,:],spm_dspec2[t,:,:]],99.99),201)
	sdlevels = np.linspace(-np.nanpercentile(np.abs(spm_dspec1[t,:,:]-spm_dspec2[t,:,:]),99),np.nanpercentile(np.abs(spm_dspec1[t,:,:]-spm_dspec2[t,:,:]),99),101)
	#
	ndspec1=np.zeros((spm_freq.shape[0],ndire.shape[0]),'f')
	ndspec1[:,1:-1]=spm_dspec1[t,:,:]
	for i in range(0,spm_freq.shape[0]):
		ndspec1[i,-1]=float((ndspec1[i,-2]+ndspec1[i,1])/2.)
		ndspec1[i,0]=float((ndspec1[i,-2]+ndspec1[i,1])/2.)

	axs[1,0].remove()
	axs[1,0] = fig.add_subplot(3, 4, 5, projection='polar')
	axs[1,0].set_theta_zero_location('N')
	axs[1,0].set_theta_direction(-1)
	axs[1,0].set_rlabel_position(-135)
	axs[1,0].set_rticks([0.1,0.15,0.20]); axs[1,0].set_rmax(1/lper)
	im = axs[1,0].contourf(theta, r, ndspec1[0:indf,:].T,slevels,cmap=plt.cm.gist_stern_r,norm=colors.PowerNorm(gamma=0.5), extend="max")
	axs[1,0].set_title('WW3Spec '+stname+', '+trun1+' '+pd.to_datetime(wtime[:][t]).strftime('%Y/%m/%d %H')+'Z',size=11) 
	del im
	# 
	ndspec2=np.zeros((spm_freq.shape[0],ndire.shape[0]),'f')
	ndspec2[:,1:-1]=spm_dspec2[t,:,:]
	for i in range(0,spm_freq.shape[0]):
		ndspec2[i,-1]=float((ndspec2[i,-2]+ndspec2[i,1])/2.)
		ndspec2[i,0]=float((ndspec2[i,-2]+ndspec2[i,1])/2.)

	axs[1,1].remove()
	axs[1,1] = fig.add_subplot(3, 4, 6, projection='polar')
	axs[1,1].set_theta_zero_location('N')
	axs[1,1].set_theta_direction(-1)
	axs[1,1].set_rlabel_position(-135)
	axs[1,1].set_rticks([0.1,0.15,0.20]); axs[1,1].set_rmax(0.27)
	im = axs[1,1].contourf(theta, r, ndspec2[0:indf,:].T,slevels,cmap=plt.cm.gist_stern_r,norm=colors.PowerNorm(gamma=0.5), extend="max")
	axs[1,1].set_title('WW3Spec '+stname+', '+trun2+' '+pd.to_datetime(wtime[:][t]).strftime('%Y/%m/%d %H')+'Z',size=11) 
	cax = axs[1,1].inset_axes([1.13, 0.2, 0.05, 0.6], transform=axs[1,1].transAxes)
	cbar = plt.colorbar(im, ax=axs[1,1], cax=cax)
	tick_locator = ticker.MaxNLocator(nbins=6); cbar.locator = tick_locator; cbar.update_ticks()
	del im
	#
	axs[1,2].remove()
	axs[1,2] = fig.add_subplot(3, 4, 7, projection='polar')
	axs[1,2].set_theta_zero_location('N')
	axs[1,2].set_theta_direction(-1)
	axs[1,2].set_rlabel_position(-135)
	axs[1,2].set_rticks([0.1,0.15,0.20]); axs[1,2].set_rmax(1/lper)
	im = axs[1,2].contourf(theta, r, ndspec2[0:indf,:].T - ndspec1[0:indf,:].T ,sdlevels,cmap=plt.cm.RdBu_r,extend="both")
	axs[1,2].set_title('WW3Spec '+stname+', '+trun1+' - '+trun2,size=10) 
	cax = axs[1,2].inset_axes([1.13, 0.2, 0.05, 0.6], transform=axs[1,2].transAxes)
	cbar = plt.colorbar(im, ax=axs[1,2], cax=cax)
	tick_locator = ticker.MaxNLocator(nbins=6); cbar.locator = tick_locator; cbar.update_ticks()
	del im
	# -----------------
	# Wave Spectra Comparison, WW3 with Buoy -----------------
	axs[1,3].remove()
	axs[1,3] = fig.add_subplot(3, 4, 8, projection='rectilinear')
	axs[1,3].scatter(spm_freq[0:indf],spm_dmspec1[t,0:indf], color='b', marker='x', label=trun1, zorder=2)
	axs[1,3].scatter(spm_freq[0:indf],spm_dmspec2[t,0:indf], color='r', marker='+', label=trun2, zorder=2)
	axs[1,3].grid(linewidth=0.5, color='grey', alpha=0.5, linestyle='--', zorder=1)
	axs[1,3].set_xlabel('Frequency (Hz)', fontsize=12); axs[1,3].set_ylabel('Direction (°)', fontsize=12) 
	axs[1,3].set_title('DirSpec '+stname+' '+pd.to_datetime(wtime[:][t]).strftime('%Y/%m/%d %H')+'Z')
	indt=np.where(spm_time[t]==spo_time)
	if size(indt)>0:
		axs[1,3].scatter(spo_freq[0:indfb],np.nanmean(spo_dmspec[indt[0][:],0:indfb],axis=0), color='k', marker='s', label='buoy', zorder=1)

	axs[1,3].legend(loc='best')

	indt=np.where(spm_time[t]==spo_time)
	# *******   Spectral Interpolation Block  *******
	# Spectral Interpolation, Buoy and WW3, taking into account the spectral/frequency bands
	# Buoy
	spi_freq = np.arange(np.max([np.nanmin(spo_freq),np.nanmin(spm_freq)]), np.min([np.nanmax(spo_freq),np.nanmax(spm_freq)]), 0.0001) # very high resolution to guarantee the integral gives exactly the same origional Hs
	spoi_pspec = np.copy(spi_freq)*0.
	for i in range(1,spi_freq.shape[0]):		
		idf=np.where( ((spi_freq[i]-spo_freq1)>=0) & ((spi_freq[i]-spo_freq2)<=0) )[0][0]
		spoi_pspec[i] = 100.*np.nanmean(spo_pspec[indt[0][:],:],axis=0)[idf]*(spi_freq[i]-spi_freq[i-1])/spo_dfreq[idf]
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

	if size(indt)>0:
		axs[2,3].plot(spi_freq[1:3650],spoi_pspec[1:3650], color='grey', linestyle='-',linewidth=0.5, label='buoy', zorder=2)
		axs[2,3].legend(loc='best')
		axs[2,3].fill_between(spi_freq[1:3650], 0.,spoi_pspec[1:3650], color='silver', alpha=0.7,label='buoy', zorder=1)

	axs[2,3].grid(linewidth=0.5, color='grey', alpha=0.5, linestyle='--', zorder=1)
	axs[2,3].set_xlabel('Frequency (Hz)', fontsize=12); axs[2,3].set_ylabel('Power Spectrum (m$^2$/Hz)', fontsize=12) 
	axs[2,3].set_ylim(ymin = -0.001)
	axs[2,3].set_title('PowerSpec '+stname+' '+pd.to_datetime(wtime[:][t]).strftime('%Y/%m/%d %H')+'Z')
	axs[2,3].axis('tight')

	# -----------------
	# Time-Series -----------------
	indtst=np.where(tsm_time==wtime[t])[0][0]
	axs[2,0].remove(); axs[2,1].remove(); axs[2,2].remove()
	iaux=np.intersect1d(tsm_time[indtst-60:indtst+61], np.array(tso_time, dtype='datetime64[h]'), assume_unique=False, return_indices=True)

	axs[2,0] = fig.add_subplot(3, 4, 9, projection='rectilinear')
	axs[2,0].plot_date(tsm_time[indtst-60:indtst+61],tsm_hs1[indtst-60:indtst+61],color='b', linestyle='--',marker='',linewidth=2.0, label=trun1, zorder=3)
	axs[2,0].plot_date(tsm_time[indtst-60:indtst+61],tsm_hs2[indtst-60:indtst+61],color='r', linestyle='-.',marker='',linewidth=2.0, label=trun2, zorder=3)
	if size(iaux)>0:
		axs[2,0].plot_date(tso_time[np.min(iaux[2]):np.max(iaux[2])+1],tso_hs[np.min(iaux[2]):np.max(iaux[2])+1],'k.',label='buoy', zorder=2)
	axs[2,0].xaxis.set_major_formatter( DateFormatter('%b%d') ); axs[2,0].fmt_xdata = DateFormatter('%b%d')
	axs[2,0].axvline(x=wtime[t],color='grey', zorder=1)
	axs[2,0].grid(linewidth=0.5, color='grey', alpha=0.5, linestyle='--', zorder=1)
	axs[2,0].legend(loc='best', fontsize=9)
	axs[2,0].set_xlabel('Date', fontsize=12); axs[2,0].set_ylabel('hs (m)', fontsize=12) 
	axs[2,0].axis('tight')
	axs[2,0].set_xlim( tsm_time[indtst-60:indtst+61].min(), tsm_time[indtst-60:indtst+61].max() )
	#
	axs[2,1] = fig.add_subplot(3, 4, 10, projection='rectilinear')
	axs[2,1].plot_date(tsm_time[indtst-60:indtst+61],tsm_tm1[indtst-60:indtst+61],color='b', linestyle='--',marker='',linewidth=2.0, label=trun1, zorder=3)
	axs[2,1].plot_date(tsm_time[indtst-60:indtst+61],tsm_tm2[indtst-60:indtst+61],color='r', linestyle='-.',marker='',linewidth=2.0, label=trun2, zorder=3)
	if size(iaux)>0:
		axs[2,1].plot_date(tso_time[np.min(iaux[2]):np.max(iaux[2])+1],tso_tm[np.min(iaux[2]):np.max(iaux[2])+1],'k.',label='buoy', zorder=2)
	axs[2,1].xaxis.set_major_formatter( DateFormatter('%b%d') ); axs[2,1].fmt_xdata = DateFormatter('%b%d')
	axs[2,1].axvline(x=wtime[t],color='grey', zorder=1)
	axs[2,1].grid(linewidth=0.5, color='grey', alpha=0.5, linestyle='--', zorder=1)
	axs[2,1].legend(loc='best', fontsize=9)
	axs[2,1].set_xlabel('Date', fontsize=12); axs[2,1].set_ylabel('tm (s)', fontsize=12) 
	axs[2,1].axis('tight')
	axs[2,1].set_xlim( tsm_time[indtst-60:indtst+61].min(), tsm_time[indtst-60:indtst+61].max() )
	#
	axs[2,2] = fig.add_subplot(3, 4, 11, projection='rectilinear')
	axs[2,2].plot_date(tsm_time[indtst-60:indtst+61],tsm_dm1[indtst-60:indtst+61],color='b', linestyle='--',marker='',linewidth=2.0, label=trun1, zorder=3)
	axs[2,2].plot_date(tsm_time[indtst-60:indtst+61],tsm_dm2[indtst-60:indtst+61],color='r', linestyle='-.',marker='',linewidth=2.0, label=trun2, zorder=3)
	if size(iaux)>0:
		axs[2,2].plot_date(tso_time[np.min(iaux[2]):np.max(iaux[2])+1],tso_dm[np.min(iaux[2]):np.max(iaux[2])+1],'k.',label='buoy', zorder=2)
	axs[2,2].xaxis.set_major_formatter( DateFormatter('%b%d') ); axs[2,2].fmt_xdata = DateFormatter('%b%d')
	axs[2,2].axvline(x=wtime[t],color='grey', zorder=1)
	axs[2,2].grid(linewidth=0.5, color='grey', alpha=0.5, linestyle='--', zorder=1)
	axs[2,2].legend(loc='best', fontsize=9)
	axs[2,2].set_xlabel('Date', fontsize=12); axs[2,2].set_ylabel('dm (°)', fontsize=12) 
	axs[2,2].axis('tight')
	axs[2,2].set_xlim( tsm_time[indtst-60:indtst+61].min(), tsm_time[indtst-60:indtst+61].max() )

	fig.canvas.draw() # https://github.com/SciTools/cartopy/issues/1207
	fig.tight_layout()
	plt.savefig('ComparisonBoard_'+stname+'_'+np.str(pd.to_datetime(wtime[:][t]).strftime('%Y%m%d%H'))+'_'+trun1+'_'+trun2+'.png', dpi=300, facecolor='w', edgecolor='w',
		orientation='portrait', papertype=None, format='png',transparent=False, pad_inches=0.1)
		
	plt.close('all'); del axs, fig, indtst, iaux

