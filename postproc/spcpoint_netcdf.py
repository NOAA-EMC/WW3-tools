# This program reads WAVEWATCHIII point output netcdf file _spec.nc from ww3_ounp
# It organizes and creates figures of power spectra (1D) and directional spectra (2D)
# It additionally divides each directional spectrum (2D) in five frequency band to be further ploted using PLEDS/WW3 (Parente (1999); a time evolution plot of frequency bands)
# Customize and select your options in lines 39 to 70.
# v1.2, 17/02/2016
#

# Python Libraries. Pay attention to the pre-requisites and libraries
import matplotlib
matplotlib.use('Agg')
import os
import pylab
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
import time
from time import strptime, gmtime, strftime
from calendar import timegm
# Palette and colors for plotting the figures
from mpl_toolkits.basemap import cm
colormap = cm.GMT_polar
palette = plt.cm.gist_stern_r
# palette.set_bad('aqua', 10.0)
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
from matplotlib.mlab import *
import sys

# Decision of ploting output figures and(or) writing output files. If equal to 0 the output will not be created. If equal to 1 the output will be created
# Power Spectra (1D) Plot 
pspdecp=0
# Directional Wave spectra (2D) Plot
dspdecp=0
# Directional Wave spectra (2D) Plot plus Power spectra and info
cspdecp=1
# PREP PLEDS output file
ppoutf=1
#
# lowest period (upper limit frequency) for the directional wave spectra (2D) polat plot
lper=3.5
# lowest period (upper limit frequency) for the power spectra (2D) plot
lper2=1.0
# Prep PLEDS (spectral partitioning) ---------------------------------------------------
# Index limits for each frequency band. Take a look at the file directions_frequencies.txt and choose the frequency interval ans its indexes.
# 1) initialFreq-la, 2) la-lb, 3) lb-lc, 4) lc-ld, 5) ld-finalFreq  : Total of 5 bands 
la=9
lb=16
lc=23
ld=30
# Output path
lpath=str(sys.argv[1])
sysid=str(sys.argv[2])
datestr=str(sys.argv[3])
# Path and Name of spectra output file
# /home//forecastop/output/ww3_OP.20160214/POINTOUTPUT/ww3.20160214_spec.nc
#aux=strftime("%Y-%m-%d %H:%M:%S", gmtime()); 
#datestr=aux[0:4]+aux[5:7]+aux[8:10]  # year month day
fname=lpath+'/'+sysid+'.'+datestr+'/netcdf/ww3.'+datestr[0:8]+'_spec.nc'
fbname='./bnames.txt'
# final path for saving figures
fpath=lpath+'/'+sysid+'.'+datestr+'/plots/points/'
tpath=lpath+'/'+sysid+'.'+datestr+'/tables/'
# --------------------------------------------------------------------------

otf=0; c=1 # Wait up to 18 hours for the output file fname
while otf==0 and c<226:
	if os.path.isfile(fname):
		otf=1
	else:
		c=c+1
		if c>10:
			time.sleep(300)


if otf==1:

	try:
		# open output file, Check if it is properly opened
		nfile = netCDF4.Dataset(fname)
	except:
		print('Problems trying to open file '+fname)
	else:
		# Read in buoy names
		with open(fbname) as f:
			bname = f.readlines()

		# efth(time, station, frequency, direction)
		dspec=nfile.variables['efth'][:,:,:,:]

		# number of time outputs 
		nt=dspec.shape[0]
		# number of point outputs 
		npo=dspec.shape[1]
		# number of directions
		nd=dspec.shape[3]
		# number of frequencies
		nf=dspec.shape[2]

		# directions
		dire=nfile.variables['direction'][:]
		# frequencies
		freq=nfile.variables['frequency'][:]
		# wind intensity and wind direction
		wnds=nfile.variables['wnd'][:,:]
		wndd=nfile.variables['wnddir'][:,:]
		# days since 1990-01-01T00:00:00Z to seconds since 1979-01-01T00:00:00Z
		atime=nfile.variables['time'][:]
		ntime=np.zeros((nt),'d')
		for t in range(0,nt):
		    ntime[t]=atime[t]*3600*24 + timegm( strptime('19900101','%Y%m%d') )

		# water depth (constant in time)
		depth1=nfile.variables['dpt'][:,:]
		depth=np.zeros((depth1.shape[1]),'f')
		for i in range(0,depth1.shape[1]):
		    depth[i]=depth1[:,i].mean()

		del depth1
		stationid=nfile.variables['station'][:]
		stationname=nfile.variables['station_name'][:,:]
		lon=nfile.variables['longitude'][0,:]
		lat=nfile.variables['latitude'][0,:]

		nfile.close()
		# -----------------------------------------------------


		# organizing directions  -----
		adspec=np.copy(dspec); inddire=int(find(dire==min(dire)))
		for t in range(0,nt):
		    for p in range(0,npo):
		        adspec[t,p,:,0:nd-(inddire+1)]=dspec[t,p,:,(inddire+1):nd]
		        adspec[t,p,:,nd-(inddire+1):nd]=dspec[t,p,:,0:(inddire+1)]
		        for i in range(0,nd):
		            dspec[t,p,:,i]=adspec[t,p,:,nd-i-1]
			    
		        adspec[t,p,:,0:int(nd/2)]=dspec[t,p,:,int(nd/2):nd]
		        adspec[t,p,:,int(nd/2):nd]=dspec[t,p,:,0:int(nd/2)]
		        dspec[t,p,:,:]=adspec[t,p,:,:]

		dire=np.sort(dire)
		# -----
		# for the 2D polar plot:
		indf=int(find(abs(freq-(1/lper))==min(abs(freq-(1/lper)))))
		indf2=int(find(abs(freq-(1/lper2))==min(abs(freq-(1/lper2)))))
		ndire=np.zeros((dire.shape[0]+2),'f'); ndire[1:-1]=dire[:]; ndire[0]=0; ndire[-1]=360
		angle = np.radians(ndire)
		r, theta = np.meshgrid(freq[0:indf], angle);
		# ---


		sdate=np.zeros((nt,4),'i')
		for t in range(0,nt):
		    sdate[t,0]=int(time.gmtime(ntime[t])[0])
		    sdate[t,1]=int(time.gmtime(ntime[t])[1])
		    sdate[t,2]=int(time.gmtime(ntime[t])[2])
		    sdate[t,3]=int(time.gmtime(ntime[t]+10)[3])

		# DF in frequency (dfim) . Nelson Violante and Fred Ostritz       
		fretab=np.zeros((nf),'f');  dfim=np.zeros((nf),'f')	
		fre1 = freq[0] ; fretab[0] = fre1
		co = freq[(nf-1)]/freq[(nf-2)] 
		dfim[0] = (co-1) * np.pi / nd * fretab[0] 
		for ifre in range(1,nf-1):
			fretab[ifre] = fretab[ifre-1] * co 
			dfim[ifre] = (co-1) * np.pi / nd * (fretab[ifre]+fretab[ifre-1]) 

		fretab[nf-1] = fretab[nf-2]*co
		dfim[nf-1] = (co-1) * np.pi / nd * fretab[(nf-2)]       
		# ------------------ 

		# intialization of variables that will be used to calculate the energy, Hs and Dp of each band(PLEDS).
		ef=np.zeros((nt,npo,5),'f') 
		dpf=np.zeros((nt,npo,5),'f')     
		pws=np.zeros((nt,npo,nf),'f')
		pwst=np.zeros((nt,npo,nf),'f') 
		hs=np.zeros((nt,npo,5),'f') 
		hstot=np.zeros((nt,npo),'f') 
		indd=np.zeros((nt,npo,5),'i') 
		# -----

		# Main loop in time
		for t in range(0,nt):
			# second loop: point outputs 
			for p in range(0,npo):

				# PREP PLEDS (Campos&Parente 2009) -------------------------------------------------
				for il in range(0,nf):	
		     			# efth(time, station, frequency, direction)
					pwst[t,p,il]=sum(dspec[t,p,il,:])

				pwst[t,p,:]=pwst[t,p,:]*dfim[:]	
		
				# Total energy at each frequency band
				ef[t,p,0]=sum(pwst[t,p,0:la])
				ef[t,p,1]=sum(pwst[t,p,la:lb])
				ef[t,p,2]=sum(pwst[t,p,lb:lc])
				ef[t,p,3]=sum(pwst[t,p,lc:ld])
				ef[t,p,4]=sum(pwst[t,p,ld:])
				# Significant wave height at each frequency band
				hs[t,p,0]=4.01*np.sqrt(ef[t,p,0])
				hs[t,p,1]=4.01*np.sqrt(ef[t,p,1])
				hs[t,p,2]=4.01*np.sqrt(ef[t,p,2])
				hs[t,p,3]=4.01*np.sqrt(ef[t,p,3])
				hs[t,p,4]=4.01*np.sqrt(ef[t,p,4])
		
				# Just for checking the total significant wave height. It must be equal to 4.01*np.sqrt(sum(pwst[t,p,:]))
				hstot[t,p]=np.sqrt(hs[t,p,0]**2 + hs[t,p,1]**2 + hs[t,p,2]**2 + hs[t,p,3]**2 + hs[t,p,4]**2)
				
				# Directional index of the maximum energy at each band				
				indd[t,p,0] = int(dspec[t,p,0:la,:].argmax()/dspec[t,p,0:la,:].shape[0])
				indd[t,p,1] = int(dspec[t,p,la:lb,:].argmax()/dspec[t,p,la:lb,:].shape[0]) 			
				indd[t,p,2] = int(dspec[t,p,lb:lc,:].argmax()/dspec[t,p,lb:lc,:].shape[0]) 
				indd[t,p,3] = int(dspec[t,p,lc:ld,:].argmax()/dspec[t,p,lc:ld,:].shape[0]) 
				indd[t,p,4] = int(dspec[t,p,ld:nf,:].argmax()/dspec[t,p,ld:nf,:].shape[0]) 	
				
				dpf[t,p,0] = dire[indd[t,p,0]]
				dpf[t,p,1] = dire[indd[t,p,1]]
				dpf[t,p,2] = dire[indd[t,p,2]]
				dpf[t,p,3] = dire[indd[t,p,3]]
				dpf[t,p,4] = dire[indd[t,p,4]]	
				# --------------------------------------------------------------------------------------------


				# SPECTRAL PLOTS 

				# Power Spectra Plot -------------------------------------------------
				if pspdecp==1:
					fig=plt.figure(figsize=(8,6))
					plt.plot(freq[0:indf2],pwst[t,p,0:indf2],'k',linewidth=3.0)
					zpwst = np.zeros((pwst[t,p,0:indf2].shape[0]),'f')
					plt.fill_between(freq[0:indf2],pwst[t,p,0:indf2], where=pwst[t,p,0:indf2]>=zpwst, interpolate=True, color='0.5')
					plt.grid()
					plt.xlabel('Frequency (Hz)')
					plt.ylabel('Power Spectrum (m^2/Hz)')         
					title('WW3 Run '+datestr+'], Wave Power Spectrum '+str(sdate[t,0]).zfill(4)+'/'+str(sdate[t,1]).zfill(2)+'/'+str(sdate[t,2]).zfill(2)+' '+str(sdate[t,3]).zfill(2)+'Z', fontsize=10)		
					plt.figtext(0.63,0.88,"Water Depth: "+repr(depth[p])[0:6]+" m",fontsize=9)
					plt.figtext(0.63,0.85,"Wind Speed: "+repr(wnds[t,p])[0:4]+" m/s",fontsize=9)
					plt.figtext(0.63,0.82,"Wind Direc: "+repr(wndd[t,p])[0:4]+" degrees",fontsize=9)
					plt.figtext(0.63,0.79,"Significant Wave Height: "+repr(4.01*np.sqrt(sum(pwst[t,p,:])))[0:4]+" m",fontsize=9)
					plt.savefig(fpath+'pspec_p'+bname[p][:-1]+'_t'+str(t).zfill(3)+'.jpg', dpi=None, facecolor='w', edgecolor='w',
					orientation='portrait', papertype=None, format='jpg',
					transparent=False, bbox_inches=None, pad_inches=0.1)
					plt.close()
				# -------------------------------------------------------------------


				# Directional Spectra 2D Plot ---------------------------------
				if dspdecp==1:

					# Polar Plot ===========================================
					ndspec=np.zeros((freq.shape[0],ndire.shape[0]),'f')
					ndspec[:,1:-1]=dspec[t,p,:,:]
					for i in range(0,freq.shape[0]):
						ndspec[i,-1]=(dspec[t,p,i,0]+dspec[t,p,i,-1])/2
						ndspec[i,0]=(dspec[t,p,i,0]+dspec[t,p,i,-1])/2

					fig=plt.figure(figsize=(8,6))
					gs = gridspec.GridSpec(10,10)
					sp = plt.subplot(gs[1:,:], projection='polar')
					sp.set_theta_zero_location('N')
					sp.set_theta_direction(-1)
					if dspec[t,p,:,:].max()<6.0:
						levels=[0.02,0.04,0.08,0.1,0.2,0.3,0.5,0.8,1,2,3,4,5,6]
						plt.contourf(theta, r, ndspec[0:indf,:].T,levels,cmap=palette,norm=colors.Normalize(vmin=0.02,vmax=6.0,clip=False), extent=[-3,3,-3,3])
					else:
						levels=np.linspace(0.1,dspec[t,p,:,:].max(), num=12, endpoint=True, retstep=False)
						plt.contourf(theta, r, ndspec[0:indf,:].T,levels.round(1),cmap=palette,norm=colors.Normalize(vmin=0.1,vmax=dspec[t,p,:,:].max(),clip=False), extent=[-3,3,-3,3])

					plt.figtext(0.15,0.93,'WW3 [Run '+datestr+'], Directional Wave Spectrum at '+bname[p][:-1]+' '+str(sdate[t,0]).zfill(4)+'/'+str(sdate[t,1]).zfill(2)+'/'+str(sdate[t,2]).zfill(2)+' '+str(sdate[t,3]).zfill(2)+'Z', fontsize=13)
					ax = plt.gca()
					pos = ax.get_position()
					l, b, w, h = pos.bounds
					cax = plt.axes([l+w+0.01, b+0.01, 0.03, h-0.01]) # setup colorbar axes.
					plt.colorbar(cax=cax) # draw colorbar
					plt.axes(ax)  # make the original axes current again
					plt.savefig(fpath+'dspec_'+bname[p][:-1]+'_t'+str(t).zfill(3)+'.jpg', dpi=None, facecolor='w', edgecolor='w',
					orientation='portrait', papertype=None, format='jpg',
					transparent=False, bbox_inches=None, pad_inches=0.1)
					plt.close()
					# ===================================================


				# Directional Spectra 2D Plot plus Power spectra and info  ---------------------------
				if cspdecp==1:

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

					plt.figtext(0.2,0.93,'Wave Spectrum at '+namep+',  '+str(sdate[t,0]).zfill(4)+'/'+str(sdate[t,1]).zfill(2)+'/'+str(sdate[t,2]).zfill(2)+' '+str(sdate[t,3]).zfill(2)+'Z', fontsize=14)
					ax = plt.gca()
					ax.set_rmin(-0.04)
					pos = ax.get_position()
					l, b, w, h = pos.bounds
					cax = plt.axes([l+w+0.040, b+0.01, 0.02, h-0.01]) # setup colorbar axes.
					cbar=plt.colorbar(cax=cax, orientation='vertical',format='%3.1f');cbar.ax.tick_params(labelsize=8)
					plt.colorbar(cax=cax) # draw colorbar
					cbar.ax.set_ylabel('Normalized by max[E(f,Dir)] = '+str(round(np.amax(dspec[t,:,:]),4))+'m2/Hz/Deg', rotation=270, fontsize=9, labelpad=10)
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

					plt.figtext(0.05,0.84,"Latitude: "+str(round(abs(lat),2))+"$^\circ$"+fllat+"    Longitude: "+str(round(abs(lon),2))+"$^\circ$"+fllon,fontsize=12)
					plt.figtext(0.05,0.75,"Water Depth: "+str(round(depth,1))+" m",fontsize=13)
					plt.figtext(0.05,0.68,"Significant Wave Height: "+str(round(4.01*np.sqrt(sum(pwst[t,:])),2))+" m",fontsize=13)
					plt.figtext(0.05,0.61,"Wind Speed: "+str(round(wnds[t],2))+" m/s",fontsize=13)
					plt.figtext(0.05,0.54,"Wind Direc: "+str(round(wndd[t]))+"$^\circ$",fontsize=13)
					plt.savefig(fpath+'/dspec_full_'+namep+'_t'+str(t).zfill(3)+'.jpg', dpi=None, facecolor='w', edgecolor='w',
					orientation='portrait', papertype=None, format='jpg',
					transparent=False, bbox_inches=None, pad_inches=0.1)
					plt.close()
					# ===================================================

				# -----------------------------------------

				plt.close('all')

		# Writing output text files
		# PREP PLEDS output file ---------------------------------------------
		if ppoutf==1:
			for p in range(0,npo):

				vf=file(tpath+'prepPLEDS_'+bname[p][:-1]+'.txt','w')
				vf.write('% WW3 [Run '+datestr+'], Input file for PLEDS-WW3 plots (Campos&Parente 2009)')
				vf.write('\n')
				vf.write('% Position:  Lat '+repr(lat[p])[0:6]+'   Lon '+repr(lon[p])[0:6])
				vf.write('\n')
				vf.write('% depth:  '+repr(depth[p])+' m')
				vf.write('\n')
				vf.write('% Freq Band (s):    1) '+repr(1/freq[0])[0:5]+'-'+repr(1/freq[la])[0:5]+'        2) '+repr(1/freq[la])[0:5]+'-'+repr(1/freq[lb])[0:5]+'        3) '+repr(1/freq[lb])[0:5]+'-'+repr(1/freq[lc])[0:5]+'        4) '+repr(1/freq[lc])[0:5]+'-'+repr(1/freq[ld])[0:5]+'        5) '+repr(1/freq[ld])[0:5]+'-'+repr(1/freq[nf-1])[0:5])         
				vf.write('\n')
				vf.write('% Date(Y M D H)     En(1) Hs(1) Dp(1)     En(2) Hs(2)  Dp(2)    En(3) Hs(3) Dp(3)     En(4) Hs(4) Dp(4)     En(5) Hs(5) Dp(5)')
				vf.write('\n')
				np.savetxt(vf,np.transpose([sdate[:,0],sdate[:,1],sdate[:,2],sdate[:,3],ef[:,p,0],hs[:,p,0],dpf[:,p,0],ef[:,p,1],hs[:,p,1],dpf[:,p,1],ef[:,p,2],hs[:,p,2],dpf[:,p,2],ef[:,p,3],hs[:,p,3],dpf[:,p,3],ef[:,p,4],hs[:,p,4],dpf[:,p,4]]),fmt="%4i %2i %2i %2i %12.5f %4.2f %5.1f %10.5f %4.2f %5.1f %10.5f %4.2f %5.1f %10.5f %4.2f %5.1f %10.5f %4.2f %5.1f",delimiter='/t') 
				vf.close()

		# -------------------------------------------------------------------

