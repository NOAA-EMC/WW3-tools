#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ww3fields.py

VERSION AND LAST UPDATE:
 v1.0  04/04/2022
 v1.1  04/27/2022

PURPOSE:
 Wave Field Map plots of WAVEWATCHIII results using cartopy. 
 It reads both netcdf and grib2 format.

USAGE:
 Mandatory Inputs: fileName and VariableName
 Examples (from linux/terminal command line):
  python3 ww3fields.py ww3gefs.20160921_field.nc hs
  python3 ww3fields.py ww3gefs.20160921_field.grib2 swh
  python3 ww3fields.py ww3gefs.20160921_field.nc t01 2 
  python3 ww3fields.py ww3gefs.20160921_field.nc phs1 1 [-10,65] [-150,10]
  nohup python3 ww3fields.py ww3gefs.20160921_field.nc hs 1 [10,60] [-100,10] >> nohup_ww3fields_20160921_hs.out 2>&1 &
 Additional optional Inputs after the Mandatory Inputs: 
  skipTime (plot at lower time resolution), 1(plot everything), 2(skip one time) etc
  min/max latitude (select an specific smaller domain)
  min/max longitude (select an specific smaller domain)

OUTPUT:
 png figures of the wave fields of selected variable of choice.
 If you want to change the resolution (for publications), edit savefig

DEPENDENCIES:
 See dependencies.py and the imports below.

AUTHOR and DATE:
 04/04/2022: Ricardo M. Campos, first version.
 04/27/2022: Ricardo M. Campos, Ali Abdolali, Matthew Masarik, Saeideh Banihashemi:
   new grids added, including tripolar and unstructured grids. New global and polar projections

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import matplotlib
matplotlib.use('Agg')
import xarray as xr
import netCDF4 as nc
import numpy as np
from pylab import *
from calendar import timegm
from time import strptime
from datetime import datetime, timezone
import matplotlib.pyplot as plt
import sys
import pandas as pd
import cartopy
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
from matplotlib import ticker
# import pickle
import warnings; warnings.filterwarnings("ignore")

palette = plt.cm.jet

# Input arguments
sk=1; slat=0; slon=0
if len(sys.argv) < 3 :
	sys.exit(' Two inputs must be provided: fileName and VariableName')
elif len(sys.argv) == 3:
	fname=np.str(sys.argv[1]); wvar=np.str.lower(sys.argv[2])
elif len(sys.argv) == 4:
	fname=np.str(sys.argv[1]); wvar=np.str.lower(sys.argv[2]); sk=np.int(sys.argv[3])
elif len(sys.argv) == 5:
	fname=np.str(sys.argv[1]); wvar=np.str.lower(sys.argv[2]); sk=np.int(sys.argv[3])
	slat=sys.argv[4]; slat=slat[1:-1]; slat=np.array([np.float(np.str(slat).split(',')[0]),np.float(np.str(slat).split(',')[1])])
elif len(sys.argv) == 6:
	fname=np.str(sys.argv[1]); wvar=np.str.lower(sys.argv[2]); sk=np.int(sys.argv[3])
	slat=sys.argv[4]; slat=slat[1:-1]; slat=np.array([np.float(np.str(slat).split(',')[0]),np.float(np.str(slat).split(',')[1])])
	slon=sys.argv[5]; slon=slon[1:-1]; slon=np.array([np.float(np.str(slon).split(',')[0]),np.float(np.str(slon).split(',')[1])])
elif len(sys.argv) > 6:
	sys.exit(' Too many inputs')

# ----- READ DATA ----- 
if np.str(fname).split('.')[-1] == 'grib2' or np.str(fname).split('.')[-1] == 'grb2': 
	# grib2 format
	ds = xr.open_dataset(fname, engine='cfgrib')
	wtime = np.array(ds.time.values + ds.step.values )
	if wvar=='hs':
		wvar=np.str('swh')

	if (wvar in list(ds.keys()))==False:
		sys.exit(' Variable name not included in the file. You can use ncdump -h filename to see variable names.')

	if size(ds[wvar].shape)==3:
		# Structured
		gstr=2
		wdata = np.array(np.flip(ds[wvar].values[:,:,:],1))
	elif size(ds[wvar].shape)==2:
		# Unstructured
		gstr=1
		wdata = np.array(ds[wvar].values)
	else:
		sys.exit(' Unexpected file shape.')

	units_wdata = np.str(ds[wvar].units)
	lat = np.array(ds.latitude.values); lon = np.array(ds.longitude.values)
	if gstr==2 and size(lat.shape)==1:
		lat = np.sort(lat)

	ds.close(); del ds
	# -----------

else:
	# netcdf format
	f=nc.Dataset(fname)
	if (wvar in list(f.variables.keys()))==False:
		sys.exit(' Variable name not included in the file. You can use ncdump -h filename to see variable names.')

	wdata = np.array(f.variables[wvar])
	units_wdata = np.str(f.variables[wvar].units)
	lat = np.array(f.variables['latitude']); lon = np.array(f.variables['longitude'])
	if size(wdata.shape)==3 and size(lat.shape)==1:
		# Structured
		gstr=2
	elif size(wdata.shape)==2 or size(lat.shape)==2:
		# Unstructured
		gstr=1

	# time
	auxt = np.array(f.variables['time'][:]*24*3600 + timegm( strptime(np.str(f.variables['time'].units).split(' ')[2][0:4]+'01010000', '%Y%m%d%H%M') )).astype('double')
	f.close(); del f
	
	wtime=[]
	for i in range(0,auxt.shape[0]):
		wtime = np.append(wtime,datetime.fromtimestamp(auxt[i], timezone.utc))

	del auxt

	wdata[wdata>360]=np.nan

if size(wdata.shape)==3 and size(lat.shape)==2:
	lon[lon>180]=lon[lon>180]-360.

if np.any(slat):
	slat=np.array(np.sort(slat))
else:
	slat=np.array([np.nanmin(lat),np.nanmax(lat)])

if np.any(slon):
	slon=np.array(np.sort(slon))
else:
	slon=np.array([np.nanmin(lon),np.nanmax(lon)])


# cbar levels
extdm=1
if "dir" in wvar or "dp" in wvar or "wlv" in wvar or "u" in wvar or "v" in wvar:
	levels = np.linspace(np.nanmin(wdata),np.nanmax(wdata),101)
	extdm=0
elif "fp" in wvar:
	levels = np.linspace(0.,np.nanpercentile(wdata,99.995),101)
else:
	levels = np.linspace(np.nanmin(wdata),np.nanpercentile(wdata,99.995),101)

print("ww3fields.py maps, file: "+fname+",  field: "+wvar)
# loop time
for t in range(wtime[::sk].shape[0]):

	# Robinson projection if global grid
	if np.diff(slat)>150 and np.diff(slon)>350 and gstr==2:
		# ax=plt.axes(projection=ccrs.Mollweide())
		ax=plt.axes(projection=ccrs.Robinson())
		gl = ax.gridlines(draw_labels=True,x_inline=False, y_inline=False, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
		gl.ylocator = ticker.MultipleLocator(15)
		gl.xlabel_style = {'size': 7, 'color': 'black', 'rotation': 0, 'rotation_mode': 'anchor'}
		gl.ylabel_style = {'size': 7, 'color': 'black', 'rotation': 0, 'rotation_mode': 'anchor'}
		data=wdata[::sk,:,:][t,:,:]
		adata, alon = add_cyclic_point(data, coord=lon)
		alat=lat 

	# Polar projection
	elif slat.max()>87 and slat.min()>30:
		ax=plt.axes(projection=ccrs.NorthPolarStereo())
		ax.set_extent([-180, 180, slat.min(), 90], crs=ccrs.PlateCarree())
		# ax.tricontourf(lon[ind],lat[ind],wdata[::sk,:][t,:][ind],levels,transform=ccrs.PlateCarree())
		gl = ax.gridlines(draw_labels=True,x_inline=False, y_inline=False, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
		gl.ylocator = ticker.MultipleLocator(10)
		gl.xlabel_style = {'size': 7, 'color': 'black', 'rotation': 0, 'rotation_mode': 'anchor'}
		gl.ylabel_style = {'size': 7, 'color': 'black', 'rotation': 0, 'rotation_mode': 'anchor'}
		alon=lon; alat=lat
		adata=wdata[::sk,:][t,:]

	# Miller projection
	else:
		ax = plt.axes(projection=ccrs.PlateCarree())
		ax.set_extent([slon.min(),slon.max(),slat.min(),slat.max()], crs=ccrs.PlateCarree())
		gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--', zorder=3)
		gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
		alon=lon; alat=lat
		adata=wdata[::sk,:][t,:]

	if gstr==1:
		# Unstructured
		ind=np.where(np.isnan(wdata[::sk,:][t,:])==False)
		if extdm==1:	
			if slat.max()>87 and slat.min()>30:
				cs=ax.tricontourf(lon[ind],lat[ind],wdata[::sk,:][t,:][ind],levels,cmap=palette,extend="max", zorder=1,transform=ccrs.PlateCarree())
			else:
				cs=ax.tricontourf(lon[ind],lat[ind],wdata[::sk,:][t,:][ind],levels,cmap=palette,extend="max", zorder=1)
		else:
			if slat.max()>87 and slat.min()>30:
				cs=ax.tricontourf(lon[ind],lat[ind],wdata[::sk,:][t,ind],levels,cmap=palette,zorder=1,transform=ccrs.PlateCarree())
			else:
				cs=ax.tricontourf(lon[ind],lat[ind],wdata[::sk,:][t,ind],levels,cmap=palette,zorder=1)
	else:
		# Structured
		if extdm==1:
			cs=ax.contourf(alon,alat,adata,levels,cmap=palette,extend="max", zorder=1,transform = ccrs.PlateCarree())
		else:
			cs=ax.contourf(alon,alat,adata,levels,cmap=palette,zorder=1,transform = ccrs.PlateCarree())

	ax.add_feature(cartopy.feature.OCEAN,facecolor=("white"))
	ax.add_feature(cartopy.feature.LAND,facecolor=("lightgrey"), edgecolor='grey',linewidth=0.5, zorder=2)
	ax.add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1, zorder=3)
	ax.coastlines(resolution='110m', color='grey',linewidth=0.5, linestyle='-', alpha=1, zorder=4)
	plt.title(wvar+' ('+units_wdata+')    '+pd.to_datetime(wtime[::sk][t]).strftime('%Y/%m/%d %H:%M')+'Z') 
	plt.tight_layout()
	ax = plt.gca()
	pos = ax.get_position()
	l, b, w, h = pos.bounds
	cax = plt.axes([l+0.07, b-0.07, w-0.12, 0.025]) # setup colorbar axes.
	cbar=plt.colorbar(cs,cax=cax, orientation='horizontal'); cbar.ax.tick_params(labelsize=10)
	tick_locator = ticker.MaxNLocator(nbins=7); cbar.locator = tick_locator; cbar.update_ticks()
	plt.axes(ax)  # make the original axes current again
	plt.tight_layout()
	# plt.savefig('wfields_'+wvar+'_'+np.str(pd.to_datetime(wtime[::sk][t]).strftime('%Y%m%d%H'))+'.eps', format='eps', dpi=200)
	plt.savefig('wfields_'+wvar+'_'+np.str(pd.to_datetime(wtime[::sk][t]).strftime('%Y%m%d%H%M'))+'.png', dpi=200, facecolor='w', edgecolor='w',
		orientation='portrait', papertype=None, format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
	
	# The pickle line below allow users to save the figure and edit later, using pickle. 
	#   to use this option, uncomment the line below and the initial import pickle
	# pickle.dump(ax, open('wfields_'+wvar+'_'+np.str(pd.to_datetime(wtime[::sk][t]).strftime('%Y%m%d%H'))+'.pickle', 'wb'))
	plt.close('all'); del ax



# For gif animation using .png figures:
# convert -delay 15 -loop 0 wfields_*.png wfields.gif

