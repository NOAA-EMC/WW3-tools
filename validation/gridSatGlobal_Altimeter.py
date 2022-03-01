import numpy as np
from matplotlib.mlab import *
from pylab import *
import os
import netCDF4 as nc
import pyresample
from mpl_toolkits.basemap import shiftgrid
import time
from calendar import timegm
import warnings
warnings.filterwarnings("ignore")
tstart = time.time()
fnetcdf="NETCDF4"

# nohup ipython3 gridSat_altimeter.py 0 >> nohup_sat0.txt 2>&1 &
npcs=5 # number of procs for parallelization
pia=10 # power of initial array 10**pia that will be used to allocate satellite data
dlim=25000.
maxti=1800.
dirs='/media/data/observations/satellite/altimeter/AODN_altm'
wpath='/home/dataproc/2collocation'
datelim='2021123123'
# Satellite missions selected, pick one (s)!
s=np.int(sys.argv[1]) # satellite ID for satellite mission selection.
sdname=np.array(['JASON3','JASON2','CRYOSAT2','JASON1','HY2','SARAL','SENTINEL3A','ENVISAT','ERS1','ERS2','GEOSAT','GFO','TOPEX'])
sname=np.array(['JASON-3','JASON-2','CRYOSAT-2','JASON-1','HY-2','SARAL','SENTINEL-3A','ENVISAT','ERS-1','ERS-2','GEOSAT','GFO','TOPEX'])

# Quality Control parameters
max_swh_rms = 1.5  # Max RMS of the band significant wave height
max_sig0_rms = 0.8 # Max RMS of the backscatter coefficient
max_swh_qc = 2.0 # Max SWH Ku band quality control
hsmax=20.; wspmax=60.
min_swh_numval = np.array([17,17,17,17,17,17,17,17,17,17,-inf,3,7])

# weight function
def wf(pdist):
	a=(1 - pdist / (dlim+1))
	return (abs(a)+a)/2

# Mask with lat lon arrays you want to collocate the altimeter data into.
f=nc.Dataset('gridInfo_GEFS.nc')
latm=f.variables['latitude'][:]; lonm=f.variables['longitude'][:]
mask=f.variables['mask'][:,:]; f.close(); del f
mask,lonm = shiftgrid(180.,mask,lonm,start=False)

# Sat files (squares) considering the domain (latm lonm) of interest, for the AODN file names
auxlat=np.array(np.arange(-90.,90.+1.,1)).astype('int')
auxlon=np.array(np.arange(0.,360.+1,1)).astype('int')

# Pyresample target points (latm lonm for valid water points).
NEWLON, NEWLAT = meshgrid(lonm,latm)
ind=np.where(mask>0); flatm = np.array(NEWLAT[ind]); flonm = np.array(NEWLON[ind]); del NEWLON, NEWLAT; del ind
targ_def = pyresample.geometry.SwathDefinition(lons=flonm,lats=flatm)

# Read and allocate satellite data into arrays
ast=np.double(np.zeros((10**pia),'d')); aslat=np.zeros((10**pia),'f'); aslon=np.zeros((10**pia),'f');
ahsku=np.zeros((10**pia),'f'); ahskucal=np.zeros((10**pia),'f'); ahsc=np.zeros((10**pia),'f');
awnd=np.zeros((10**pia),'f'); awndcal=np.zeros((10**pia),'f'); asig0kunstd=np.zeros((10**pia),'f');
aswhkunobs=np.zeros((10**pia),'f'); aswhkunstd=np.zeros((10**pia),'f'); aswhkuqc=np.zeros((10**pia),'f')
ii=0
for j in auxlat:
	for k in auxlon:

		if j>=0:
			hem='N'
		else:
			hem='S'

		try: 
			fu=nc.Dataset(dirs+'/'+sdname[s]+'/IMOS_SRS-Surface-Waves_MW_'+sname[s]+'_FV02_'+str(np.abs(j)).zfill(3)+hem+'-'+str(k).zfill(3)+'E-DM00.nc')
		except:
			print(dirs+'/'+sdname[s]+'/IMOS_SRS-Surface-Waves_MW_'+sname[s]+'_FV02_'+str(np.abs(j)).zfill(3)+hem+'-'+str(k).zfill(3)+'E-DM00.nc does not exist'); vai=0
		else:
			st=np.double(fu.variables['TIME'][:])
			if size(st)>10:
				slat=fu.variables['LATITUDE'][:]
				slon=fu.variables['LONGITUDE'][:]
				hsc=fu.variables['SWH_C'][:]
				wnd=fu.variables['WSPD'][:]
				wndcal=fu.variables['WSPD_CAL'][:]
				try: 
					hsku=fu.variables['SWH_KU'][:]
					hskucal=fu.variables['SWH_KU_CAL'][:]
					sig0kunstd=fu.variables['SIG0_KU_std_dev'][:]
					swhkunobs=fu.variables['SWH_KU_num_obs'][:]
					swhkunstd=fu.variables['SWH_KU_std_dev'][:]
					swhkuqc=fu.variables['SWH_KU_quality_control'][:]
				except:
					print(' error reading KU, pick KA')
					hsku=fu.variables['SWH_KA'][:]
					hskucal=fu.variables['SWH_KA_CAL'][:]
					sig0kunstd=fu.variables['SIG0_KA_std_dev'][:]
					swhkunobs=fu.variables['SWH_KA_num_obs'][:]
					swhkunstd=fu.variables['SWH_KA_std_dev'][:]
					swhkuqc=fu.variables['SWH_KA_quality_control'][:]


				if ii+size(st) <= ast.shape[0] :
					if (st.shape[0]==wnd.shape[0]) & (slat.shape[0]==slon.shape[0]) & (hsku.shape[0]==hskucal.shape[0]) :	
						ast[ii:ii+st.shape[0]]=np.array(st).astype('double')
						aslat[ii:ii+st.shape[0]]=np.array(slat).astype('float')
						aslon[ii:ii+st.shape[0]]=np.array(slon).astype('float')
						ahsku[ii:ii+st.shape[0]]=np.array(hsku).astype('float')
						ahskucal[ii:ii+st.shape[0]]=np.array(hskucal).astype('float')
						ahsc[ii:ii+st.shape[0]]=np.array(hsc).astype('float')
						awnd[ii:ii+st.shape[0]]=np.array(wnd).astype('float')
						awndcal[ii:ii+st.shape[0]]=np.array(wndcal).astype('float')
						asig0kunstd[ii:ii+st.shape[0]]=np.array(sig0kunstd).astype('float')
						aswhkunobs[ii:ii+st.shape[0]]=np.array(swhkunobs).astype('float')
						aswhkunstd[ii:ii+st.shape[0]]=np.array(swhkunstd).astype('float')
						aswhkuqc[ii:ii+st.shape[0]]=np.array(swhkuqc).astype('float')
						ii=ii+st.shape[0]

				else:
					sys.exit('Small array to allocate the satellite data! Increase the power of initial array (pia)')

				del st,slat,slon,hsku,hskucal,hsc,wnd,sig0kunstd,swhkunobs,swhkunstd,swhkuqc
				fu.close(); del fu

print(' Done reading and allocating satellite data '+sdname[s])
del ii

# Quality Control Check ----
indq = np.where( (aswhkunstd<=max_swh_rms) & (asig0kunstd<=max_sig0_rms) & (aswhkunobs>=min_swh_numval[s]) & (aswhkuqc<=max_swh_qc) & (ahsku>0.01) & (ahsku<hsmax) & (awnd>0.01) & (awnd<wspmax) & (ast<=float(timegm( time.strptime(datelim, '%Y%m%d%H') ))))     
del asig0kunstd,aswhkunobs,aswhkunstd,aswhkuqc
if np.any(indq):
	if size(indq)>10:
		ii=0
		ast=np.double(np.copy(ast[indq[0]]))
		ast=np.double(np.copy(ast)*24.*3600.+float(timegm( time.strptime('1985010100', '%Y%m%d%H') )))
		aslat=np.copy(aslat[indq[0]]); aslon=np.copy(aslon[indq[0]])
		ahsku=np.copy(ahsku[indq[0]]); ahskucal=np.copy(ahskucal[indq[0]]); ahsc=np.copy(ahsc[indq[0]])
		awnd=np.copy(awnd[indq[0]]); awndcal=np.copy(awndcal[indq[0]])
		# Collocated Arrays:
		# final hourly time array. Combined with flatm and flonm, it represents the reference for sat collocation
		atime = np.array(np.arange(np.double((np.round(ast.min()/3600.)*3600.)-3600.),np.double((np.round(ast.max()/3600.)*3600.)+3600.)+1,3600.)).astype('double')
		#
		ftime=np.double(np.zeros((ast.shape[0]*2),'d')); flat=np.zeros((ast.shape[0]*2),'f'); flon=np.zeros((ast.shape[0]*2),'f')
		fhsku=np.zeros((ast.shape[0]*2),'f'); stdhsku=np.zeros((ast.shape[0]*2),'f'); counthsku=np.zeros((ast.shape[0]*2),'f')
		fhskucal=np.zeros((ast.shape[0]*2),'f'); fhsc=np.zeros((ast.shape[0]*2),'f') 
		fwnd=np.zeros((ast.shape[0]*2),'f'); fwndcal=np.zeros((ast.shape[0]*2),'f')

		# into the regular grid with kd tree
		for t in range(0,atime.shape[0]):
			indt = np.where( abs(ast[:]-atime[t]) < maxti )
			if size(indt)>2:
				prlon=np.copy(aslon[indt[0]]); prlon[prlon>180.]=prlon[prlon>180.]-360.
				orig_def = pyresample.geometry.SwathDefinition(lons=prlon, lats=aslat[indt[0]]); del prlon
				# By distance function wf
				auxfhsku, auxstdhsku, auxcounthsku = pyresample.kd_tree.resample_custom(orig_def,ahsku[indt[0]],targ_def,radius_of_influence=dlim,weight_funcs=wf,fill_value=0,with_uncert=True,nprocs=npcs)
				auxfhskucal = pyresample.kd_tree.resample_custom(orig_def,ahskucal[indt[0]],targ_def,radius_of_influence=dlim,weight_funcs=wf,fill_value=0,nprocs=npcs)
				auxfhsc = pyresample.kd_tree.resample_custom(orig_def,ahsc[indt[0]],targ_def,radius_of_influence=dlim,weight_funcs=wf,fill_value=0,nprocs=npcs)
				auxfwnd = pyresample.kd_tree.resample_custom(orig_def,awnd[indt[0]],targ_def,radius_of_influence=dlim,weight_funcs=wf,fill_value=0,nprocs=npcs)
				auxfwndcal = pyresample.kd_tree.resample_custom(orig_def,awndcal[indt[0]],targ_def,radius_of_influence=dlim,weight_funcs=wf,fill_value=0,nprocs=npcs)
				# print('   - ok '+repr(t))
				indpqq = np.where( (auxfhskucal>0.01) & (auxfwnd>0.01) & (auxfhskucal<hsmax) & (auxfwnd<wspmax) )[0]
				# allocate data into final array
				ftime[ii:ii+size(indpqq)] = np.array(np.zeros(size(indpqq),'d')+atime[t]).astype('double')
				flat[ii:ii+size(indpqq)] = np.array(flatm[indpqq]).astype('float')
				flon[ii:ii+size(indpqq)] = np.array(flonm[indpqq]).astype('float')
				fhsku[ii:ii+size(indpqq)] = np.array(auxfhsku[indpqq]).astype('float')
				stdhsku[ii:ii+size(indpqq)] = np.array(auxstdhsku[indpqq]).astype('float')
				counthsku[ii:ii+size(indpqq)] = np.array(auxcounthsku[indpqq]).astype('float')
				fhskucal[ii:ii+size(indpqq)] = np.array(auxfhskucal[indpqq]).astype('float')
				fhsc[ii:ii+size(indpqq)] = np.array(auxfhsc[indpqq]).astype('float')
				fwnd[ii:ii+size(indpqq)] = np.array(auxfwnd[indpqq]).astype('float')
				fwndcal[ii:ii+size(indpqq)] = np.array(auxfwndcal[indpqq]).astype('float')
				ii=ii+size(indpqq)

				del auxfhsku, auxstdhsku, auxcounthsku, auxfhskucal, auxfhsc, auxfwnd, auxfwndcal, indpqq
				
			del indt
			print('PyResample kdtree, hourly time, '+repr(t))

		del ast, aslat, aslon, ahsku, ahskucal, ahsc, awnd, awndcal, indq, atime

		print(' '); print(sdname[s]+' Done')

		# Final quality control (double check)
		ind=np.where( (fhsku<0.01) | (fhsku>hsmax) )
		if np.any(ind):
			fhsku[ind[0]]=np.nan; del ind

		ind=np.where( (fhskucal<0.01) | (fhskucal>hsmax) )
		if np.any(ind):
			fhskucal[ind[0]]=np.nan; del ind

		ind=np.where( (fhsc<0.01) | (fhsc>hsmax) )
		if np.any(ind):
			fhsc[ind[0]]=np.nan; del ind

		ind=np.where( (fwnd<0.01) | (fwnd>wspmax) )
		if np.any(ind):
			fwnd[ind[0]]=np.nan; del ind

		ind=np.where( (fwndcal<0.01) | (fwndcal>wspmax) )
		if np.any(ind):
			fwndcal[ind[0]]=np.nan; del ind

		indf=np.where( (ftime>0.) & (fhsku>=0.0) )
		ftime=np.array(ftime[indf[0]]).astype('double')
		flat=np.array(flat[indf[0]]); flon=np.array(flon[indf[0]])
		fhsku=np.array(fhsku[indf[0]]); fhskucal=np.array(fhskucal[indf[0]])
		stdhsku=np.array(stdhsku[indf[0]]); counthsku=np.array(counthsku[indf[0]])
		fhsc=np.array(fhsc[indf[0]])
		fwnd=np.array(fwnd[indf[0]]); fwndcal=np.array(fwndcal[indf[0]])
		print(sdname[s]+' .  Array Size '+np.str(size(indf))); print(' ')
		del indf

		# Save netcdf
		ncfile = nc.Dataset('AltimeterGridded_'+sdname[s]+'_onGEFSgrid.nc', "w", format=fnetcdf) 
		ncfile.history="AODN Altimeter data on GEFSv12 0.25X0.25 grid." 
		# ncfile.author = "Ricardo Campos"
		# create  dimensions.
		ncfile.createDimension('time' , ftime.shape[0] )
		# create variables.
		vflat = ncfile.createVariable('latitude',np.dtype('float32').char,('time'))
		vflon = ncfile.createVariable('longitude',np.dtype('float32').char,('time'))
		vft = ncfile.createVariable('stime',np.dtype('float64').char,('time'))
		vfhsku = ncfile.createVariable('hsku',np.dtype('float32').char,('time'))
		vstdhsku = ncfile.createVariable('stdhsku',np.dtype('float32').char,('time'))
		vcounthsku = ncfile.createVariable('counthsku',np.dtype('float32').char,('time'))
		vfhskucal = ncfile.createVariable('hskucal',np.dtype('float32').char,('time'))
		vfhsc = ncfile.createVariable('hsc',np.dtype('float32').char,('time'))
		vfwnd = ncfile.createVariable('wnd',np.dtype('float32').char,('time'))
		vfwndcal = ncfile.createVariable('wndcal',np.dtype('float32').char,('time'))
		# Assign units
		vft.units = 'seconds since 1970-01-01 00:00:00'
		vflat.units = 'degrees_north' ; vflon.units = 'degrees_east'
		vfwnd.units = 'm/s' ; vfwndcal.units = 'm/s'
		vfhsku.units = 'm'; vfhskucal.units = 'm'; vfhsc.units = 'm'
		# Allocate Data
		vflat[:] = flat; vflon[:] = flon; vft[:] = ftime
		vfhsku[:] = fhsku; vstdhsku[:] = stdhsku; vcounthsku[:] = counthsku
		vfhskucal[:] = fhskucal; vfhsc[:] = fhsc
		vfwnd[:] = fwnd; vfwndcal[:] = fwndcal
		ncfile.close()
		print(' ')
		print('netcdf ok ')

		del vft, vflat, vflon, vfhsku, vstdhsku, vcounthsku, vfhskucal, vfhsc, vfwnd, vfwndcal


