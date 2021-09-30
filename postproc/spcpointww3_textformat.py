# This program reads WAVEWATCHIII point output from NCEP/NOAA operational forecasts
# check the points at https://polar.ncep.noaa.gov/waves/viewer.shtml?-gfswave-

# Python Libraries. Pay attention to the pre-requisites and libraries
import matplotlib
matplotlib.use('Agg')
import os
import pylab
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
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
import warnings
warnings.filterwarnings("ignore")

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
#
fname=str(sys.argv[1]) # input file name
nt=np.int(sys.argv[2]) # number of time outputs
fpath=str(sys.argv[3]) # final path, output figs
tpath=str(sys.argv[4]) # final path, output table
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
		fp = open(fname)
	except:
		print('Problems trying to open file '+fname)
	else:
		# Header --------------------------------------
		cabc=fp.readline()
		cabc=cabc.strip().split()
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

		nl=int(floor((nf*nd)/7.));rnl=int(round(7*((float(nf*nd)/7)-nl)))
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
		r, theta = np.meshgrid(freq[0:indf], angle);
		# ---

		sdate=np.zeros((nt,4),'i')
		for t in range(0,nt):
		    sdate[t,0]=int(time.gmtime(ntime[t])[0])
		    sdate[t,1]=int(time.gmtime(ntime[t])[1])
		    sdate[t,2]=int(time.gmtime(ntime[t])[2])
		    sdate[t,3]=int(time.gmtime(ntime[t]+10)[3])

		# DF in frequency (dfim)     
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
		ef=np.zeros((nt,5),'f') 
		dpf=np.zeros((nt,5),'f')     
		pws=np.zeros((nt,nf),'f')
		pwst=np.zeros((nt,nf),'f') 
		hs=np.zeros((nt,5),'f') 
		hstot=np.zeros((nt),'f') 
		indd=np.zeros((nt,5),'i') 
		# -----

		# Main loop in time
		for t in range(0,nt):

			# PREP PLEDS (Campos&Parente 2009) -------------------------------------------------
			for il in range(0,nf):	
	     			# efth(time, frequency, direction)
				pwst[t,il]=sum(dspec[t,il,:])

			pwst[t,:]=pwst[t,:]*dfim[:]	
		
			# Total energy at each frequency band
			ef[t,0]=sum(pwst[t,0:la])
			ef[t,1]=sum(pwst[t,la:lb])
			ef[t,2]=sum(pwst[t,lb:lc])
			ef[t,3]=sum(pwst[t,lc:ld])
			ef[t,4]=sum(pwst[t,ld:])
			# Significant wave height at each frequency band
			hs[t,0]=4.01*np.sqrt(ef[t,0])
			hs[t,1]=4.01*np.sqrt(ef[t,1])
			hs[t,2]=4.01*np.sqrt(ef[t,2])
			hs[t,3]=4.01*np.sqrt(ef[t,3])
			hs[t,4]=4.01*np.sqrt(ef[t,4])
		
			# Just for checking the total significant wave height. It must be equal to 4.01*np.sqrt(sum(pwst[t,:]))
			hstot[t]=np.sqrt(hs[t,0]**2 + hs[t,1]**2 + hs[t,2]**2 + hs[t,3]**2 + hs[t,4]**2)
				
			# Directional index of the maximum energy at each band				
			indd[t,0] = int(dspec[t,0:la,:].argmax()/dspec[t,0:la,:].shape[0])
			indd[t,1] = int(dspec[t,la:lb,:].argmax()/dspec[t,la:lb,:].shape[0]) 			
			indd[t,2] = int(dspec[t,lb:lc,:].argmax()/dspec[t,lb:lc,:].shape[0]) 
			indd[t,3] = int(dspec[t,lc:ld,:].argmax()/dspec[t,lc:ld,:].shape[0]) 
			indd[t,4] = int(dspec[t,ld:nf,:].argmax()/dspec[t,ld:nf,:].shape[0]) 	
				
			dpf[t,0] = dire[indd[t,0]]
			dpf[t,1] = dire[indd[t,1]]
			dpf[t,2] = dire[indd[t,2]]
			dpf[t,3] = dire[indd[t,3]]
			dpf[t,4] = dire[indd[t,4]]	
			# --------------------------------------------------------------------------------------------


			# SPECTRAL PLOTS 

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

			plt.close('all')

		# Writing output text files
		# PREP PLEDS output file ---------------------------------------------
		vf=open(tpath+'/prepPLEDS_'+namep+'.txt','w')
		vf.write('% WW3 spec output, Input file for PLEDS-WW3 plots')
		vf.write('\n')
		vf.write('% Position:  Lat '+repr(lat)[0:6]+'   Lon '+repr(lon)[0:6])
		vf.write('\n')
		vf.write('% depth:  '+repr(depth)+' m')
		vf.write('\n')
		vf.write('% Freq Band (s):    1) '+repr(1/freq[0])[0:5]+'-'+repr(1/freq[la])[0:5]+'        2) '+repr(1/freq[la])[0:5]+'-'+repr(1/freq[lb])[0:5]+'        3) '+repr(1/freq[lb])[0:5]+'-'+repr(1/freq[lc])[0:5]+'        4) '+repr(1/freq[lc])[0:5]+'-'+repr(1/freq[ld])[0:5]+'        5) '+repr(1/freq[ld])[0:5]+'-'+repr(1/freq[nf-1])[0:5])         
		vf.write('\n')
		vf.write('% Date(Y M D H)     En(1) Hs(1) Dp(1)     En(2) Hs(2)  Dp(2)    En(3) Hs(3) Dp(3)     En(4) Hs(4) Dp(4)     En(5) Hs(5) Dp(5)')
		vf.write('\n')
		np.savetxt(vf,np.transpose([sdate[:,0],sdate[:,1],sdate[:,2],sdate[:,3],ef[:,0],hs[:,0],dpf[:,0],ef[:,1],hs[:,1],dpf[:,1],ef[:,2],hs[:,2],dpf[:,2],ef[:,3],hs[:,3],dpf[:,3],ef[:,4],hs[:,4],dpf[:,4]]),fmt="%4i %2i %2i %2i %12.5f %4.2f %5.1f %10.5f %4.2f %5.1f %10.5f %4.2f %5.1f %10.5f %4.2f %5.1f %10.5f %4.2f %5.1f",delimiter='/t') 
		vf.close()

		# -------------------------------------------------------------------

