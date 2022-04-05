"""
ww3fields.py

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
  nohup python3 ww3fields.py ww3gefs.20160921_field.nc hs 1 [10,60] [-100,10] >> nohup_ww3fields_20160921_hs.txt 2>&1 &
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

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import matplotlib
matplotlib.use('Agg')
import xarray as xr
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import sys
import pandas as pd
import cartopy.crs as ccrs
import cartopy
from matplotlib import ticker
# import pickle
import warnings
warnings.filterwarnings("ignore")

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

	wdata = np.array(np.flip(ds[wvar].values[:,:,:],1))

else:
	# netcdf format
	ds = xr.open_dataset(fname)
	wtime = np.array(ds.time.values)
	wdata = np.array(ds[wvar].values[:,:,:])

units_wdata = np.str(ds[wvar].units)
lat = np.sort(np.array(ds.latitude.values[:]))
lon = np.array(ds.longitude.values[:])
ds.close(); del ds
# -----------
if np.any(slat):
	slat=np.sort(slat)
else:
	slat=np.array([np.nanmin(lat),np.nanmax(lat)])

if np.any(slon):
	slon=np.sort(slon)
	if slon.min()<-180. :
		sys.exit(' Longitude below -180. Keep the longitude standard: 0to360 or -180to180 degrees.')

else:
	slon=np.array([np.nanmin(lon),np.nanmax(lon)])

# cbar levels
levels = np.linspace(np.nanmin(wdata),np.nanpercentile(wdata,99.99),101)
print("ww3fields.py maps, file: "+fname+",  field: "+wvar)
# loop time
for t in range(wdata[::sk,:,:].shape[0]):
	fig, ax = plt.subplots()
	ax = plt.axes(projection=ccrs.PlateCarree()) 
	ax.set_extent([slon.min(),slon.max(),slat.min(),slat.max()], crs=ccrs.PlateCarree())  
	gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
	gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
	plt.contourf(lon,lat,wdata[::sk,:,:][t,:,:],levels,cmap=palette,extend="max", zorder=2)
	ax.add_feature(cartopy.feature.OCEAN,facecolor=("white"))
	ax.add_feature(cartopy.feature.LAND,facecolor=("lightgrey"), edgecolor='grey',linewidth=0.5)
	ax.add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1)
	ax.coastlines(resolution='110m', color='grey',linewidth=0.5, linestyle='-', alpha=1)
	plt.title(wvar+' ('+units_wdata+')    '+pd.to_datetime(wtime[::sk][t]).strftime('%Y/%m/%d %H')+'Z') 
	fig.tight_layout()
	ax = plt.gca()
	pos = ax.get_position()
	l, b, w, h = pos.bounds
	cax = plt.axes([l+0.07, b-0.07, w-0.15, 0.025]) # setup colorbar axes.
	cbar=plt.colorbar(cax=cax, orientation='horizontal'); cbar.ax.tick_params(labelsize=10)
	tick_locator = ticker.MaxNLocator(nbins=7); cbar.locator = tick_locator; cbar.update_ticks()
	plt.axes(ax)  # make the original axes current again
	fig.tight_layout()
	# plt.savefig('wfields_'+wvar+'_'+np.str(pd.to_datetime(wtime[::sk][t]).strftime('%Y%m%d%H'))+'.eps', format='eps', dpi=200)
	plt.savefig('wfields_'+wvar+'_'+np.str(pd.to_datetime(wtime[::sk][t]).strftime('%Y%m%d%H'))+'.png', dpi=200, facecolor='w', edgecolor='w',
		orientation='portrait', papertype=None, format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
	
	# pickle.dump(fig, open('wfields_'+wvar+'_'+np.str(pd.to_datetime(wtime[::sk][t]).strftime('%Y%m%d%H'))+'.pickle', 'wb'))
	plt.close('all'); del ax, fig


