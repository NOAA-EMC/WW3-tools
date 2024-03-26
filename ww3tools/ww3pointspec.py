#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ww3pointspec.py

VERSION AND LAST UPDATE:
 v1.0  04/04/2022
 v2.0  03/25/2024

PURPOSE:
 Directional Spectrum plots of WAVEWATCHIII results using polar plots.
 It reads both netcdf and text format.

USAGE:
 Mandatory Inputs: fileName and StationName (or ID, starting with 0)
 Additional optional Inputs after the Mandatory Inputs: 
   skipTime (plot at lower time resolution) where 1(plot everything), 2(skip one time) etc
 Examples (from linux/terminal command line):
  python3 ww3pointspec.py ww3gefs.20160928_spec.nc 41002
  python3 ww3pointspec.py ww3gefs.20160928_spec.nc 41002 2
  nohup python3 ww3pointspec.py ww3gefs.20160928_spec.nc 41004 >> nohup_ww3pointspec_20160928.out 2>&1 &
 A shell script, get_ww3spec_c00.sh, can be run to test this code 
  ww3pointspec.py operationally, plotting today's forecast spectra.

OUTPUT:
 png figures of the directional and power spectrum for the specific 
  point/station selected as input.
 If you want to change the resolution (for publications), edit savefig

DEPENDENCIES:
 See setup.py and the imports below.

AUTHOR and DATE:
 04/04/2022: Ricardo M. Campos, first version.
 03/25/2024: Ricardo M. Campos, using wread.spec_ww3 to read the spectra

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import matplotlib
matplotlib.use('Agg')
import xarray as xr
import numpy as np
from pylab import *
from matplotlib.mlab import *
import matplotlib.pyplot as plt
import sys
import pandas as pd
from time import strptime, strftime
from calendar import timegm
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import wread
import warnings; warnings.filterwarnings("ignore")
# Palette and colors for plotting the figures
palette = plt.cm.gist_stern_r

# lowest period (upper limit frequency) for the directional wave spectra (2D) polat plot
lper=3.5
# lowest period (upper limit frequency) for the power spectra (1D) plot
lper2=1.0

# Input arguments
sk=1
if len(sys.argv) < 3 :
	sys.exit(' Two inputs must be provided: fileName and StationName')
elif len(sys.argv) == 3:
	fname=str(sys.argv[1]); stname=str(sys.argv[2])
elif len(sys.argv) == 4:
	fname=str(sys.argv[1]); stname=str(sys.argv[2]); sk=int(sys.argv[3])
elif len(sys.argv) > 4:
	sys.exit(' Too many inputs')

# --- READ FILE ---
result = wread.spec_ww3(fname,stname)
freq = result['freq']; pwst = result['pspec']
dire = result['theta']; dspec = result['dirspec']
nt = int(result['time'].shape[0])

# --- PLOT ---
# for the 2D polar plot:
indf=int(np.where(abs(freq-(1/lper))==min(abs(freq-(1/lper))))[0][0])
indf2=int(np.where(abs(freq-(1/lper2))==min(abs(freq-(1/lper2))))[0][0])
ndire=np.zeros((dire.shape[0]+2),'f'); ndire[1:-1]=dire[:]; ndire[0]=0; ndire[-1]=360
angle = np.radians(ndire)
r, theta = np.meshgrid(freq[0:indf], angle)
# --- 
for t in range(0,nt,sk):

	# can confirm the spectrum agrees with Hs through 4.01*sqrt(np.sum(pwst[t,:]))

	# Polar Plot ===========================================
	ndspec=np.zeros((freq.shape[0],ndire.shape[0]),'f')
	ndspec[:,1:-1]=dspec[t,:,:]/np.amax(dspec[t,:,:])
	for i in range(0,freq.shape[0]):
		ndspec[i,-1]=float((ndspec[i,-2]+ndspec[i,1])/2.)
		ndspec[i,0]=float((ndspec[i,-2]+ndspec[i,1])/2.)
  
	fig=plt.figure(figsize=(8,4.8))
	gs = gridspec.GridSpec(10,10)
	sp = plt.subplot(gs[1:,4::], projection='polar')
	sp.set_theta_zero_location('N')
	sp.set_theta_direction(-1)
	levels=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
	plt.contourf(theta, r, ndspec[0:indf,:].T,levels,cmap=palette,norm=colors.Normalize(vmin=0.02,vmax=1.0,clip=False), extent=[-3,3,-3,3])

	plt.figtext(0.2,0.93,'Wave Spectrum at '+stname+',  '+pd.to_datetime(result['date'][t]).strftime('%Y/%m/%d %H')+'Z', fontsize=14)
	ax = plt.gca()
	ax.set_rmin(-0.04)
	pos = ax.get_position()
	l, b, w, h = pos.bounds
	cax = plt.axes([l+w+0.040, b+0.01, 0.02, h-0.01]) # setup colorbar axes.
	cbar=plt.colorbar(cax=cax, orientation='vertical',format='%3.1f');cbar.ax.tick_params(labelsize=8)
	plt.colorbar(cax=cax) # draw colorbar
	cbar.ax.set_ylabel('Normalized by max[E(f,Dir)] = '+str(np.round(np.amax(dspec[t,:,:]),4))+'m2/Hz/Deg', rotation=270, fontsize=9, labelpad=10)
	plt.axes(ax)  # make the original axes current again
	plt.figtext(0.55,0.02,"Wave Direction: Coming from",fontsize=9)

	# Power spectrum
	sp2 = plt.subplot(gs[6::,0:4])
	plt.plot(freq[0:indf2],pwst[t,0:indf2],'k',linewidth=3.0)
	sp2.tick_params(labelsize=7)
	zpwst = np.zeros((pwst[t,0:indf2].shape[0]),'f')
	plt.fill_between(freq[0:indf2],pwst[t,0:indf2], where=pwst[t,0:indf2]>=zpwst, interpolate=True, color='0.5')
	plt.grid()
	plt.xlabel('Frequency (Hz)',fontsize=9)
	plt.ylabel('Power Spectrum (m$^2$/Hz)',fontsize=9)         
	# table of info
	if result['latitude']>=0.0:
		fllat='N'
	else:
		fllat='S'
	if result['longitude']>=0.0:
		fllon='E'
	else:
		fllon='W'

	plt.figtext(0.05,0.84,"Latitude: "+str(np.round(abs(result['latitude']),2))+"$^\circ$"+fllat+"    Longitude: "+str(np.round(abs(result['longitude']),2))+"$^\circ$"+fllon,fontsize=12)
	plt.figtext(0.05,0.75,"Water Depth: "+str(np.round(result['depth'],1))+" m",fontsize=13)
	plt.figtext(0.05,0.68,"Significant Wave Height: "+str(np.round(4.01*np.sqrt(sum(pwst[t,:])),2))+" m",fontsize=13)
	plt.figtext(0.05,0.61,"Wind Speed: "+str(np.round(result['wind_spd'][t],2))+" m/s",fontsize=13)
	plt.figtext(0.05,0.54,"Wind Direc: "+str(np.round(result['wind_dir'][t]))+"$^\circ$",fontsize=13)

	plt.savefig('wspectrum_'+stname+'_'+str(pd.to_datetime(result['date'][t]).strftime('%Y%m%d%H'))+'.png', dpi=200, facecolor='w', edgecolor='w',orientation='portrait', 
		format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)

	plt.close(); del ax, fig
	# =================================================== 

# For gif animation using .png figures:
# convert -delay 15 -loop 0 wspectrum_*.png wspectrum.gif

