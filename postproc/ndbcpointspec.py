#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ndbcpointspec.py

VERSION AND LAST UPDATE:
 v1.0  06/03/2022

PURPOSE:
 Directional Spectrum and Power Spectrum Plots from NDBC.
 It reads netcdf format only. For example:
 https://dods.ndbc.noaa.gov/thredds/fileServer/data/swden/42002/42002w2020.nc
 See retrieve_ndbc_nc.py which downloads NDBC data for multiple buoys.

USAGE:
 Mandatory Inputs: fileName (can include the full path)
 Additional optional Input: 
   skipTime (plot at lower time resolution) where 1(plot everything), 2(skip one time) etc
 Examples (from linux/terminal command line):
  python3 ndbcpointspec.py 41049w2021.nc
  python3 ndbcpointspec.py 42002w2020.nc 6
  python3 ndbcpointspec.py /home/user/ndbcdata/41049w2021.nc 3
  nohup python3 ndbcpointspec.py 41049w2021.nc >> nohup_ndbcpointspec_41049.out 2>&1 &

OUTPUT:
 png figures with two subplots: the directional and power spectrum 
  for the specific point/station selected as input.
 If you want to change the resolution (for publications), edit savefig

DEPENDENCIES:
 See dependencies.py and the imports below.
 It requires the python function wread.py, at postproc of WW3-tools

AUTHOR and DATE:
 06/03/2022: Ricardo M. Campos, first version.

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import matplotlib
matplotlib.use('Agg') # for backend plots, not for rendering in a window
import numpy as np
from pylab import *
from matplotlib.mlab import *
import matplotlib.pyplot as plt
import sys
import pandas as pd
import cartopy.crs as ccrs
from matplotlib import ticker
import matplotlib.colors as colors
import warnings; warnings.filterwarnings("ignore")
import wread
palette = plt.cm.jet
dpalette = plt.cm.RdBu_r

# lowest period (upper limit frequency) for the directional wave spectra (2D) polat plot
lper=4.5
# Input arguments
#  skip time steps, to alleviate memory and reduce the number of consecutive plots.
sk=1
if len(sys.argv) < 2 :
	sys.exit(' At least one input must be provided: NDBC netcdf file')
elif len(sys.argv) == 2:
	bspec=np.str(sys.argv[1])
elif len(sys.argv) == 3:
	bspec=np.str(sys.argv[1]); sk=np.int(sys.argv[2])
elif len(sys.argv) > 3:
	sys.exit(' Too many inputs')

stname=np.str(bspec).split('.nc')[0][-10:-5]
figname="NDBCspectrum_"+np.str(bspec).split('.nc')[0][-10::]+"_"

# READ NDBC data
spo_time,ignore,spo_lat,spo_lon,spo_freq,spo_dfreq,spo_pspec,spo_dmspec,ignore,spo_dire,spo_dspec = wread.spec_ndbc(bspec,sk)

# Organize directions and frequencies
ndire=np.zeros((spo_dire.shape[0]+2),'f'); ndire[1:-1]=spo_dire[:]; ndire[0]=0; ndire[-1]=360
angle = np.radians(ndire)
indfb=int(np.where(abs(spo_freq-(1/lper))==min(abs(spo_freq-(1/lper))))[0][0])
r, theta = np.meshgrid(spo_freq[0:indfb], angle)
spo_freq1=np.array(spo_freq-spo_dfreq/2); spo_freq2=np.array(spo_freq+spo_dfreq/2)
spo_time = np.array(spo_time, dtype='datetime64[h]') # round to hour

print(" Read data Ok. Starting plots ...")
for t in range(0,spo_time.shape[0]):

	if np.any(np.isnan(spo_dspec[t,:,:])==False):
	
		ndspec=np.zeros((spo_freq.shape[0],ndire.shape[0]),'f')
		ndspec[:,1:-1]=spo_dspec[t,:,:]
		for i in range(0,spo_freq.shape[0]):
			ndspec[i,-1]=float((ndspec[i,-2]+ndspec[i,1])/2.)
			ndspec[i,0]=float((ndspec[i,-2]+ndspec[i,1])/2.)
	
		nslevels = np.linspace(0.001,np.nanpercentile(spo_dspec[t,:,:],99.99),201)	
		fig, axs = plt.subplots(nrows=1,ncols=2,subplot_kw={'projection': ccrs.PlateCarree()},figsize=(12,5))
		axs[0].remove()
		axs[0] = fig.add_subplot(1, 2, 1, projection='polar')
		axs[0].set_theta_zero_location('N')
		axs[0].set_theta_direction(-1)
		axs[0].set_rlabel_position(-135)
		axs[0].set_rticks([0.1,0.15,0.20]); axs[0].set_rmax(0.27)
		im = axs[0].contourf(theta, r, ndspec[0:indfb,:].T,nslevels,cmap=plt.cm.gist_stern_r,norm=colors.PowerNorm(gamma=0.5), extend="max")
		axs[0].set_title('NDBC Spectrum  '+stname+', '+pd.to_datetime(spo_time[t]).strftime('%Y/%m/%d %H')+'Z',size=11) 
		cax = axs[0].inset_axes([1.13, 0.2, 0.05, 0.6], transform=axs[0].transAxes)
		cbar = plt.colorbar(im, ax=axs[0], cax=cax)
		tick_locator = ticker.MaxNLocator(nbins=6); cbar.locator = tick_locator; cbar.update_ticks()
		del im
		plt.figtext(0.14,0.02,"Wave Direction: Coming from",fontsize=10)

		axs[1].remove()
		axs[1] = fig.add_subplot(1, 2, 2, projection='rectilinear')
		axs[1].plot(spo_freq,spo_pspec[t,:], color='grey', linestyle='-',linewidth=1.,zorder=2)
		axs[1].plot(spo_freq[0:indfb],spo_pspec[t,0:indfb], color='dimgrey', linestyle='-',linewidth=2.,zorder=2)
		axs[1].plot(spo_freq[:],spo_pspec[t,:],'k.',linewidth=0.5,zorder=2)
		axs[1].fill_between(spo_freq[:], 0.,spo_pspec[t,:], color='silver', alpha=0.7,label='buoy', zorder=1)
		axs[1].grid(linewidth=0.5, color='grey', alpha=0.5, linestyle='--', zorder=1)
		axs[1].set_xlabel('Frequency (Hz)', fontsize=12); axs[1].set_ylabel('Power Spectrum (m$^2$/Hz)', fontsize=12) 
		axs[1].set_title('NDBC Power Spectrum  '+stname+' '+pd.to_datetime(spo_time[t]).strftime('%Y/%m/%d %H')+'Z')
		axs[1].axis('tight'); axs[1].set_ylim(ymin = -0.0001)

		fig.canvas.draw() # https://github.com/SciTools/cartopy/issues/1207
		fig.tight_layout()
		plt.savefig(figname+np.str(pd.to_datetime(spo_time[t]).strftime('%Y%m%d%H'))+'.png', dpi=300, facecolor='w', edgecolor='w',
			orientation='portrait', papertype=None, format='png',transparent=False, pad_inches=0.1)

		plt.close('all'); del axs, fig
		print("    Figure Ok. "+np.str(pd.to_datetime(spo_time[t]).strftime('%Y%m%d%H')))

print(" Done.")
