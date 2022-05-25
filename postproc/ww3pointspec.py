#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ww3pointspec.py

VERSION AND LAST UPDATE:
 v1.0  04/04/2022

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
from matplotlib.mlab import *
import matplotlib.pyplot as plt
import sys
import pandas as pd
from time import strptime, strftime
from calendar import timegm
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
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
	fname=np.str(sys.argv[1]); stname=np.str(sys.argv[2])
elif len(sys.argv) == 4:
	fname=np.str(sys.argv[1]); stname=np.str(sys.argv[2]); sk=np.int(sys.argv[3])
elif len(sys.argv) > 4:
	sys.exit(' Too many inputs')

# --- READ FILE ---
if np.str(fname).split('.')[-1] == 'nc':
	# NetCDF format
	ds = xr.open_dataset(fname)
	dspec=np.array(ds['efth'].values[:,:,:,:])
	# number of time outputs 
	nt=dspec.shape[0]
	# number of point outputs 
	npo=dspec.shape[1]
	# number of directions
	nd=dspec.shape[3]
	# number of frequencies
	nf=dspec.shape[2]

	# directions
	dire=np.array(ds['direction'].values[:])
	# frequencies
	freq=np.array(ds['frequency'].values[:])
	# DF in frequency (dfreq)
	dfreq=np.array(ds['frequency2'].values[:] - ds['frequency1'].values[:])
	# wind intensity and wind direction
	wnds=np.array(ds['wnd'].values[:,:])
	wndd=np.array(ds['wnddir'].values[:,:])
	# Time datetime64 array
	ntime=np.array(ds['time'].values[:])
	# water depth (constant in time)
	depth=np.nanmean(ds['dpt'].values[:,:],axis=0)

	auxstationname=ds['station_name'].values[:,:]; stationname=[]
	for i in range(0,auxstationname.shape[0]):
		stationname=np.append(stationname,"".join(np.array(auxstationname[i,:]).astype('str')))

	lon=np.array(np.nanmean(ds['longitude'].values[:,:],axis=0))
	lat=np.array(np.nanmean(ds['latitude'].values[:,:],axis=0))	
	ds.close(); del ds, auxstationname

	if str.isnumeric(stname) == True:
		if np.int(stname) < 1000:
			inds=np.int(stname); stname=np.str(stationname[inds])
		else:
			inds=np.where(stationname[:]==stname)
			if size(inds)>0:
				inds=np.int(inds[0][0]); stname=np.str(stationname[inds])
			else:
				sys.exit(' Station '+stname+' not included in the output file, or wrong station ID')

	else:
		inds=np.where(stationname[:]==stname)
		if size(inds)>0:
			inds=np.int(inds[0][0]); stname=np.str(stationname[inds])
		else:
			sys.exit(' Station '+stname+' not included in the output file, or wrong station ID')

	dspec=np.array(dspec[:,inds,:,:]); wnds=np.array(wnds[:,inds]); wndd=np.array(wndd[:,inds])
	lat=np.array(lat[inds]); lon=np.array(lon[inds]); depth=np.array(depth[inds])
	del inds, stationname

else:

	# Text format (only one point allowed here, same as ww3/NCEP operational)
	fp = open(fname); nt = fp.read().count(stname); fp.close(); del fp
	if nt>=1:
		# Open file and read the first parameters
		fp = open(fname)
		cabc=fp.readline(); cabc=cabc.strip().split()
		nf=int(cabc[3]) # number of frequencies
		nd=int(cabc[4]) # number of directions
		npo=int(cabc[5]) # number of point outputs 

		freq=zeros(nf,'f');dire=zeros(nd,'f')
		dspec=zeros((nt,nf,nd),'f')
		adire=zeros(dire.shape)
		adspec=zeros(dspec.shape)
		ntime=np.zeros((nt),'d')

		# Frequencies --------------------
		ncf=np.int(np.floor(nf/8));rncf=np.int(np.round(8*((float(nf)/8)-ncf)))
		k=0
		for i in range(0,ncf):
			line=fp.readline()
			line=line.strip().split()
			for j in range(0,8):
				freq[k]=float(line[j])
				k=k+1

		if rncf>0:
			line=fp.readline()
			line=line.strip().split()
			for i in range(0,rncf):
				freq[k]=float(line[i])
				k=k+1	

		# DF in frequency (dfreq)
		dfreq=np.zeros(freq.shape[0],'f')
		for i in range(0,freq.shape[0]):
			if i==0 or i==(freq.shape[0]-1):
				dfreq[i] = freq[i]*(1+ ( ((freq[-1]/freq[-2])-1)/2 )) - freq[i]
			else:
				dfreq[i] = freq[i]*(freq[-1]/freq[-2]) - freq[i]

		# Directions ---------------------
		ncd=np.int(np.floor(nd/7));rncd=np.int(np.round(7*((float(nd)/7)-ncd)))
		k=0
		for i in range(0,ncd):
			line=fp.readline()
			line=line.strip().split()
			for j in range(0,7):
				dire[k]=float(line[j])*180/pi
				k=k+1

		if rncd>0:
			line=fp.readline()
			line=line.strip().split()
			for i in range(0,rncd):
				dire[k]=float(line[i])*180/pi
				k=k+1

		nl=int(floor((nf*nd)/7.)); rnl=int(np.round(7*((float(nf*nd)/7)-nl)))
		auxs=np.zeros((nf*nd),'f')
		wnds=np.zeros((nt),'f');wndd=np.zeros((nt),'f')
		for t in range(0,nt):
				
			cabc=fp.readline(); cabc.strip().split()[0]
			ntime[t] = np.double(timegm( strptime(cabc.strip().split()[0]+cabc.strip().split()[1][0:2], '%Y%m%d%H') ))
			cabc=fp.readline(); cabc=cabc.strip().split()
			if t==0:
				namep=cabc[0][1:]
				lat=float(cabc[2]);lon=float(cabc[3])
				depth=float(cabc[4])

			wnds[t]=float(cabc[5]);wndd[t]=float(cabc[6])

			k=0
			for i in range(0,nl):
				line=fp.readline()
				line=line.strip().split()
				for j in range(0,7):
					auxs[k]=float(line[j])
					k=k+1

			if rncd>0:
				line=fp.readline()
				line=line.strip().split()
				for i in range(0,rnl):
					auxs[k]=float(line[i])
					k=k+1

			for ic in range(0,nf):
				for il in range(0,nd):
				    dspec[t,ic,il]=auxs[il*nf+ic]
					 
		fp.close(); del fp

	else:
		sys.exit(' Station '+stname+' not included in the output file')

# -----------------------------------------------------

# organizing directions  -----
adspec=np.copy(dspec); inddire=int(np.where(dire==min(dire))[0][0])
for t in range(0,nt):
        adspec[t,:,0:nd-(inddire+1)]=dspec[t,:,(inddire+1):nd]
        adspec[t,:,nd-(inddire+1):nd]=dspec[t,:,0:(inddire+1)]
        for i in range(0,nd):
            dspec[t,:,i]=adspec[t,:,nd-i-1]
	    
        adspec[t,:,0:int(nd/2)]=dspec[t,:,int(nd/2):nd]
        adspec[t,:,int(nd/2):nd]=dspec[t,:,0:int(nd/2)]
        dspec[t,:,:]=adspec[t,:,:]

dire=np.sort(dire)
# -----
# for the 2D polar plot:
indf=int(np.where(abs(freq-(1/lper))==min(abs(freq-(1/lper))))[0][0])
indf2=int(np.where(abs(freq-(1/lper2))==min(abs(freq-(1/lper2))))[0][0])
ndire=np.zeros((dire.shape[0]+2),'f'); ndire[1:-1]=dire[:]; ndire[0]=0; ndire[-1]=360
angle = np.radians(ndire)
r, theta = np.meshgrid(freq[0:indf], angle)
# --- 

# --- PLOT ---
pwst=np.zeros((nt,nf),'f')
for t in range(0,nt,sk):

	for il in range(0,nf):	
		pwst[t,il]=sum(dspec[t,il,:]*(2*np.pi)/nd)

	pwst[t,:]=pwst[t,:]*dfreq[:]
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

	plt.figtext(0.2,0.93,'Wave Spectrum at '+stname+',  '+pd.to_datetime(ntime[t]).strftime('%Y/%m/%d %H')+'Z', fontsize=14)
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
	if lat>=0.0:
		fllat='N'
	else:
		fllat='S'
	if lon>=0.0:
		fllon='E'
	else:
		fllon='W'

	plt.figtext(0.05,0.84,"Latitude: "+str(np.round(abs(lat),2))+"$^\circ$"+fllat+"    Longitude: "+str(np.round(abs(lon),2))+"$^\circ$"+fllon,fontsize=12)
	plt.figtext(0.05,0.75,"Water Depth: "+str(np.round(depth,1))+" m",fontsize=13)
	plt.figtext(0.05,0.68,"Significant Wave Height: "+str(np.round(4.01*np.sqrt(sum(pwst[t,:])),2))+" m",fontsize=13)
	plt.figtext(0.05,0.61,"Wind Speed: "+str(np.round(wnds[t],2))+" m/s",fontsize=13)
	plt.figtext(0.05,0.54,"Wind Direc: "+str(np.round(wndd[t]))+"$^\circ$",fontsize=13)

	plt.savefig('wspectrum_'+stname+'_'+np.str(pd.to_datetime(ntime[t]).strftime('%Y%m%d%H'))+'.png', dpi=200, facecolor='w', edgecolor='w',orientation='portrait', 
		papertype=None, format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)

	plt.close(); del ax, fig
	# =================================================== 

