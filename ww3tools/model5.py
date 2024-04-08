import warnings
warnings.filterwarnings("ignore")
import numpy as np
from matplotlib.mlab import *
from pylab import *
import xarray as xr
import netCDF4 as nc
import time
from time import strptime
from calendar import timegm
import wread
# netcdf format
fnetcdf="NETCDF4"

# Paths
# ndbcp="/data/buoys/NDBC/wparam"
ndbcp="/scratch2/NCEPDEV/marine/Matthew.Masarik/dat/buoys/NDBC/ncformat/wparam"
# Copernicus buoys
# copernp="/data/buoys/Copernicus/wtimeseries"
copernp="/work/noaa/marine/ricardo.campos/data/buoys/Copernicus/wtimeseries"
print('  ')


# Options of including grid and cyclone information
gridinfo=int(0); cyclonemap=int(0); wlist=[]; ftag=''; forecastds=0



if len(sys.argv) >= 4:
	gridinfo=str(sys.argv[3])
	print(' Using gridInfo '+gridinfo)




# READ DATA
print(" ")
if gridinfo!=0:
	# Grid Information
	gridmask = wread.mask(gridinfo)
	mlat=gridmask['latitude']; mlon=gridmask['longitude']
	mask=gridmask['mask']; distcoast=gridmask['distcoast']; depth=gridmask['depth']
	oni=gridmask['GlobalOceansSeas']; ocnames=gridmask['names_GlobalOceansSeas']
	hsmz=gridmask['HighSeasMarineZones']; hsmznames=gridmask['names_HighSeasMarineZones']		
	print("  GridInfo Ok. "+gridinfo)

	# Cyclone Information
	if cyclonemap!=0:
		cycloneinfo = wread.cyclonemap(cyclonemap)
		clat=cycloneinfo['latitude']; clon=cycloneinfo['longitude']
		cmap=cycloneinfo['cmap']; ctime=cycloneinfo['time']
		cinfo=np.array(cycloneinfo['info'].split(':')[1].split(';'))
		if np.array_equal(clat,mlat)==True & np.array_equal(clon,mlon)==True: 
			print("  CycloneMap Ok. "+cyclonemap)
		else:
			sys.exit(' Error: Cyclone grid and Mask grid are different.')










from wread import spec_ww3

# Check if the correct number of arguments is provided
if len(sys.argv) != 3:
    print("Usage: python run_spec_ww3.py file_names.txt station_names.txt")
    sys.exit(1)

# Read file names from the first argument
file_names = []
with open(sys.argv[1], 'r') as file:
    for line in file:
        file_names.append(line.strip())

# Read station names from the second argument
station_names = []
with open(sys.argv[2], 'r') as file:
    for line in file:
        station_names.append(line.strip())

# Call the spec_ww3 function
results = spec_ww3(file_names, station_names)

# Additional processing
mfcycle = None
stname = None
mtime = None
mhs = None
mtp = None
mwn = None
mwd = None
mfreq = None
mtm=None
mdm=None


for i, result in enumerate(results):
    at = result['time']
    if i == 0:
        mfcycle = np.array(np.zeros((at.shape[0]), 'd') + at[0]).astype('double')
        stname = np.atleast_1d(np.array(result['station_name']))
        mtime = np.copy(at)
        mhs = np.copy([result['Hs']])
        mtp = np.copy([result['Tp']])
        if 'wind_spd' in result.keys():
            mwn = np.copy([result['wind_spd']])
        else:
            mwn = np.copy(mhs) * np.nan

        if 'wind_dir' in result.keys():
            mwd = np.copy([result['wind_dir']])
        else:
            mwd = np.copy(mhs) * np.nan

        if 'freq' in result.keys():
            mfreq = np.copy([result['freq']])
        else:
            mfreq = np.copy(mhs) * np.nan

        if 'tm' in result.keys():
            mtm = np.append(mtm,[result['tm']],axis=0)
        else:
            mtm = np.copy(mhs) * np.nan

        if 'dm' in result.keys():
           mdm=np.copy([result['dm']])
        else:
           mdm=np.copy(mhs)*np.nan

        if 'dp' in result.keys():
           mdp=np.copy([result['dp']])
        else:
           mdp=np.copy(mhs)*np.nan

    else:
        if mhs.shape[1] == result['Hs'].shape[0]:
            stname = np.append(stname, np.atleast_1d(np.array(result['station_name'])))
            mtime = np.append(mtime, at)
            mhs = np.append(mhs, [result['Hs']], axis=0)
            mtp = np.append(mtp, [result['Tp']], axis=0)
            if 'wind_spd' in result.keys():
                mwn = np.append(mwn, [result['wind_spd']], axis=0)
            else:
                mwn = np.append(mwn, [np.copy(result['Hs']) * np.nan], axis=0)

            if 'wind_dir' in result.keys():
                mwd = np.append(mwd, [result['wind_dir']], axis=0)
            else:
                mwd = np.append(mwd, [np.copy(result['Hs']) * np.nan], axis=0)

            if 'freq' in result.keys():
                mfreq = np.append(mfreq, [result['freq']], axis=0)
            else:
                mfreq = np.append(mfreq, [np.copy(result['Hs']) * np.nan], axis=0)

            if 'tm' in result.keys():
                mtm=np.append(mtm,[result['tm']],axis=0)
            else:
                mtm=np.append(mtm,[np.copy(result['Hs'])*np.nan],axis=0)

            if 'dm' in result.keys():
                mdm=np.append(mdm,[result['dm']],axis=0)
            else:
                mdm=np.append(mdm,[np.copy(result['Hs'])*np.nan],axis=0)   

            if 'dp' in result.keys():
                mdp=np.append(mdp,[result['dp']],axis=0)
            else:
                mdp=np.append(mdp,[np.copy(result['Hs'])*np.nan],axis=0)

        else:
            print("   Stations in " + wlist[i] + " do not match the other tar files. Skipped " + wlist[i])

# Save proc
print('HSS:', mhs)
print('WIND:', mwn)


# BUOYS ------------------
bwind = np.zeros((np.size(stname), np.size(mtime)), 'f') * np.nan
bhs = np.zeros((np.size(stname), np.size(mtime)), 'f') * np.nan
btm = np.zeros((np.size(stname), np.size(mtime)), 'f') * np.nan
btp = np.zeros((np.size(stname), np.size(mtime)), 'f') * np.nan
bdm = np.zeros((np.size(stname), np.size(mtime)), 'f') * np.nan
bdp = np.zeros((np.size(stname), np.size(mtime)), 'f') * np.nan
lat = np.zeros(np.size(stname), 'f') * np.nan
lon = np.zeros(np.size(stname), 'f') * np.nan
# help reading NDBC buoys, divided by year
yrange = np.array(np.arange(time.gmtime(mtime.min())[0], time.gmtime(mtime.min())[0] + 1, 1)).astype('int')
# loop buoys
for b in range(0, np.size(stname)):

    ahs = []
    try:
        awm = []
        ahs = []
        atm = []
        atp = []
        adm = []
        atime = []
        for y in yrange:

            f = nc.Dataset(ndbcp + "/" + stname[b] + "h" + repr(y) + ".nc")
            if 'wave_height' in f.variables.keys():
                ahs = np.append(ahs, f.variables['wave_height'][:, 0, 0])
            elif 'hs' in f.variables.keys():
                ahs = np.append(ahs, f.variables['hs'][:, 0, 0])
            elif 'swh' in f.variables.keys():
                ahs = np.append(ahs, f.variables['swh'][:, 0, 0])

            if 'wind_spd' in f.variables.keys():
                awm = np.append(awm, f.variables['wind_spd'][:, 0, 0])
            else:
                awm = np.array(np.copy(ahs * nan))

            if 'average_wpd' in f.variables.keys():
                atm = np.append(atm, f.variables['average_wpd'][:, 0, 0])
            else:
                atm = np.array(np.copy(ahs * nan))

            if 'dominant_wpd' in f.variables.keys():
                atp = np.append(atp, f.variables['dominant_wpd'][:, 0, 0])
            else:
                atp = np.array(np.copy(ahs * nan))

            if 'mean_wave_dir' in f.variables.keys():
                adm = np.append(adm, f.variables['mean_wave_dir'][:, 0, 0])
            else:
                adm = np.array(np.copy(ahs * nan))

            if 'latitude' in f.variables.keys():
                lat[b] = f.variables['latitude'][:]
            elif 'LATITUDE' in f.variables.keys():
                lat[b] = f.variables['LATITUDE'][:]
            else:
                lat[b] = nan

            if 'longitude' in f.variables.keys():
                lon[b] = f.variables['longitude'][:]
            elif 'LONGITUDE' in f.variables.keys():
                lon[b] = f.variables['LONGITUDE'][:]
            else:
                lon[b] = nan

            atime = np.append(atime, np.array(f.variables['time'][:]).astype('double'))

            f.close()
            del f

        adp = adm * np.nan  # no peak direction available in this format

        if np.size(ahs) > 0:

            # First layer of simple quality-control
            indq = np.where((ahs > 30.) | (ahs < 0.0))
            if np.size(indq) > 0:
                ahs[indq] = np.nan
                del indq

            indq = np.where((atm > 40.) | (atm < 0.0))
            if np.size(indq) > 0:
                atm[indq] = np.nan
                del indq

            indq = np.where((atp > 40.) | (atp < 0.0))
            if np.size(indq) > 0:
                atp[indq] = np.nan
                del indq

            indq = np.where((adm > 360.) | (adm < -180.))
            if np.size(indq) > 0:
                adm[indq] = np.nan
                del indq

            indq = np.where((adp > 360.) | (adp < -180.))
            if np.size(indq) > 0:
                adp[indq] = np.nan
                del indq

            indq = np.where((awm > 50.) | (awm < 0.0))
            if np.size(indq) > 0:
                awm[indq] = np.nan
                del indq

            c = 0
            for t in range(0, np.size(mtime)):
                indt = np.where(np.abs(atime - mtime[t]) < 1800.)
                if np.size(indt) > 0:
                    if np.any(ahs[indt[0]].mask == False):
                        bhs[b, t] = np.nanmean(ahs[indt[0]][ahs[indt[0]].mask == False])
                        c = c + 1
                    if np.any(atm[indt[0]].mask == False):
                        btm[b, t] = np.nanmean(atm[indt[0]][atm[indt[0]].mask == False])
                    if np.any(atp[indt[0]].mask == False):
                        btp[b, t] = np.nanmean(atp[indt[0]][atp[indt[0]].mask == False])
                    if np.any(adm[indt[0]].mask == False):
                        bdm[b, t] = np.nanmean(adm[indt[0]][adm[indt[0]].mask == False])
                    if np.any(adp[indt[0]].mask == False):
                        bdp[b, t] = np.nanmean(adp[indt[0]][adp[indt[0]].mask == False])
                    if np.any(awm[indt[0]].mask == False):
                        bwind[b, t] = np.nanmean(awm[indt[0]][awm[indt[0]].mask == False])

                    del indt

            # print("counted "+repr(c)+" at "+stname[b])

        print("   station " + stname[b] + "  ok")
#        del ahs
    except Exception as e:
        print("Error occurred while processing station", stname[b])
        print(e)

print('bwind:',bwind)
print('bhs:',bhs)


print('  ')
# Simple quality-control (range)
ind=np.where((bhs>30.)|(bhs<0.0))
if np.size(ind)>0:
	bhs[ind]=np.nan; del ind

ind=np.where((btm>40.)|(btm<0.0))
if np.size(ind)>0:
	btm[ind]=np.nan; del ind

ind=np.where((btp>40.)|(btp<0.0))
if np.size(ind)>0:
	btp[ind]=np.nan; del ind

ind=np.where((bdm>360.)|(bdm<-180.))
if np.size(ind)>0:
	bdm[ind]=np.nan; del ind

ind=np.where((bdp>360.)|(bdp<-180.))
if np.size(ind)>0:
	bdp[ind]=np.nan; del ind

ind=np.where((bwind>50.0)|(bwind<0.0))
if np.size(ind)>0:
        bwind[ind]=np.nan; del ind

ind=np.where((mhs>30.)|(mhs<0.0))
if np.size(ind)>0:
	mhs[ind]=np.nan; del ind

ind=np.where((mtm>40.)|(mtm<0.0))
if np.size(ind)>0:
	mtm[ind]=np.nan; del ind
	
ind=np.where((mtp>40.)|(mtp<0.0))
if np.size(ind)>0:
	mtp[ind]=np.nan; del ind

ind=np.where((mdm>360.)|(mdm<-180.))
if np.size(ind)>0:
	mdm[ind]=np.nan; del ind

ind=np.where((mdp>360.)|(mdp<-180.))
if np.size(ind)>0:
	mdp[ind]=np.nan; del ind

ind=np.where((mwn>50.)|(mwn<0.0))
if np.size(ind)>0:
        mwn[ind]=np.nan; del ind

# Clean data excluding some stations. Select matchups only when model and buoy are available.
ind=np.where( (np.isnan(lat)==False) & (np.isnan(lon)==False) & (np.isnan(np.nanmean(mhs,axis=1))==False) & (np.isnan(np.nanmean(bhs,axis=1))==False) )
if np.size(ind)>0:
	stname=np.array(stname[ind[0]])
	lat=np.array(lat[ind[0]])
	lon=np.array(lon[ind[0]])
	mhs=np.array(mhs[ind[0],:])
	mtm=np.array(mtm[ind[0],:])
	mtp=np.array(mtp[ind[0],:])
	mdm=np.array(mdm[ind[0],:])
	mdp=np.array(mdp[ind[0],:])
	mwn=np.array(mwn[ind[0],:])
	bhs=np.array(bhs[ind[0],:])
	btm=np.array(btm[ind[0],:])
	btp=np.array(btp[ind[0],:])
	bdm=np.array(bdm[ind[0],:])
	bdp=np.array(bdp[ind[0],:])
	bwind=np.array(bwind[ind[0],:])
else:
	sys.exit(' Error: No matchups Model/Buoy available.')


print(" Matchups model/buoy complete. Total of "+repr(np.size(ind))+" stations/buoys avaliable."); del ind

# Processing grid and/or cyclone information
if gridinfo!=0:
	print(" Adding extra information ... ")
	alon=np.copy(lon); alon[alon<0]=alon[alon<0]+360.
	indgplat=[]; indgplon=[]
	for i in range(0,lat.shape[0]):
		# indexes nearest point.
		indgplat = np.append(indgplat,np.where( abs(mlat-lat[i])==abs(mlat-lat[i]).min() )[0][0])
		indgplon = np.append(indgplon,np.where( abs(mlon-alon[i])==abs(mlon-alon[i]).min() )[0][0])

	indgplat=np.array(indgplat).astype('int'); indgplon=np.array(indgplon).astype('int')
	pdistcoast=np.zeros(lat.shape[0],'f')*np.nan
	pdepth=np.zeros(lat.shape[0],'f')*np.nan
	poni=np.zeros(lat.shape[0],'f')*np.nan
	phsmz=np.zeros(lat.shape[0],'f')*np.nan
	for i in range(0,lat.shape[0]):
		pdistcoast[i]=distcoast[indgplat[i],indgplon[i]]
		pdepth[i]=depth[indgplat[i],indgplon[i]]
		poni[i]=oni[indgplat[i],indgplon[i]]
		phsmz[i]=hsmz[indgplat[i],indgplon[i]]

	print(" Grid Information Included.")

	# Excluding shallow water points too close to the coast (mask information not accurate)
	ind=np.where( (np.isnan(pdistcoast)==False) & (np.isnan(pdepth)==False) )
	if np.size(ind)>0:
		stname=np.array(stname[ind[0]])
		lat=np.array(lat[ind[0]])
		lon=np.array(lon[ind[0]])
		mhs=np.array(mhs[ind[0],:])
		mtm=np.array(mtm[ind[0],:])
		mtp=np.array(mtp[ind[0],:])
		mdm=np.array(mdm[ind[0],:])
		mdp=np.array(mdp[ind[0],:])
		mwn=np.array(mwn[ind[0],:])
		bhs=np.array(bhs[ind[0],:])
		btm=np.array(btm[ind[0],:])
		btp=np.array(btp[ind[0],:])
		bdm=np.array(bdm[ind[0],:])
		bdp=np.array(bdp[ind[0],:])
		bwind=np.array(bwind[ind[0],:])
		pdistcoast=np.array(pdistcoast[ind[0]])
		pdepth=np.array(pdepth[ind[0]])
		poni=np.array(poni[ind[0]])
		phsmz=np.array(phsmz[ind[0]])
	else:
		sys.exit(' Error: No matchups Model/Buoy available after using grid mask.')

	del ind

	if cyclonemap!=0:
		fcmap=np.zeros((lat.shape[0],mtime.shape[0]),'f')*np.nan
		for t in range(0,np.size(mtime)):
			# search for cyclone time index and cyclone map
			indt=np.where(np.abs(ctime-mtime[t])<5400.)
			if np.size(indt)>0:
				for i in range(0,lat.shape[0]):
					fcmap[i,t] = np.array(cmap[indt[0][0],indgplat[i],indgplon[i]])

				del indt
			else:
				print('     - No cyclone information for this time step: '+repr(t))

			# print(' Done cyclone analysis at step: '+repr(t))

		ind=np.where(fcmap<0)
		if np.size(ind)>0:
			fcmap[ind]=np.nan

		print(" Cyclone Information Included.")

# Edit format if this is forecast model data. Reshape and allocate


if forecastds>0:
	unt=np.unique(mfcycle); mxsz=1
	for i in range(0,unt.shape[0]):
		ind=np.where(mfcycle==unt[i])[0]
		mxsz=np.max([mxsz,np.size(ind)])

	for i in range(0,unt.shape[0]):
		ind=np.where(mfcycle==unt[i])[0]
		if i==0:
			nmhs=np.zeros((mhs.shape[0],unt.shape[0],mxsz),'f')*np.nan
			nmtm=np.zeros((mhs.shape[0],unt.shape[0],mxsz),'f')*np.nan
			nmtp=np.zeros((mhs.shape[0],unt.shape[0],mxsz),'f')*np.nan
			nmdm=np.zeros((mhs.shape[0],unt.shape[0],mxsz),'f')*np.nan
			nmdp=np.zeros((mhs.shape[0],unt.shape[0],mxsz),'f')*np.nan
			nmwn=np.zeros((mhs.shape[0],unt.shape[0],mxsz),'f')*np.nan
			nbhs=np.zeros((mhs.shape[0],unt.shape[0],mxsz),'f')*np.nan
			nbtm=np.zeros((mhs.shape[0],unt.shape[0],mxsz),'f')*np.nan
			nbtp=np.zeros((mhs.shape[0],unt.shape[0],mxsz),'f')*np.nan
			nbdm=np.zeros((mhs.shape[0],unt.shape[0],mxsz),'f')*np.nan
			nbdp=np.zeros((mhs.shape[0],unt.shape[0],mxsz),'f')*np.nan
			nbwind=np.zeros((mhs.shape[0],unt.shape[0],mxsz),'f')*np.nan
			nmtime=np.zeros((unt.shape[0],mxsz),'double')*np.nan
			if cyclonemap!=0:
				nfcmap=np.zeros((mhs.shape[0],unt.shape[0],mxsz),'f')*np.nan
		
		nmtime[i,0:np.size(ind)]=np.array(mtime[ind]).astype('double')
		nmhs[:,i,:][:,0:np.size(ind)]=np.array(mhs[:,ind])
		nmtm[:,i,:][:,0:np.size(ind)]=np.array(mtm[:,ind])
		nmtp[:,i,:][:,0:np.size(ind)]=np.array(mtp[:,ind])
		nmdm[:,i,:][:,0:np.size(ind)]=np.array(mdm[:,ind])
		nmdp[:,i,:][:,0:np.size(ind)]=np.array(mdp[:,ind])
		nmwn[:,i,:][:,0:np.size(ind)]=np.array(mwn[:,ind])
		nbhs[:,i,:][:,0:np.size(ind)]=np.array(bhs[:,ind])
		nbtm[:,i,:][:,0:np.size(ind)]=np.array(btm[:,ind])
		nbtp[:,i,:][:,0:np.size(ind)]=np.array(btp[:,ind])
		nbdm[:,i,:][:,0:np.size(ind)]=np.array(bdm[:,ind])
		nbdp[:,i,:][:,0:np.size(ind)]=np.array(bdp[:,ind])
		nwind[:,i,:][:,0:np.size(ind)]=np.array(bwind[:,ind])

		if cyclonemap!=0:				
			nfcmap[:,i,:][:,0:np.size(ind)]=np.array(fcmap[:,ind])


	ind=np.where( (nmhs>0.0) & (nbhs>0.0) )

else:

    larger_shape = max(mhs.shape, bhs.shape)
    padded_mhs = np.zeros(larger_shape)
    padded_bhs = np.zeros(larger_shape)
    padded_mhs[:mhs.shape[0], :mhs.shape[1]] = mhs
    padded_bhs[:bhs.shape[0], :bhs.shape[1]] = bhs


    # Create padded arrays for other variables
    padded_mtm = np.zeros(larger_shape)  # Assuming mtm has the same shape as mhs
    padded_mtp = np.zeros(larger_shape)  # Assuming mtp has the same shape as mhs
    padded_mdm = np.zeros(larger_shape)  # Assuming mdm has the same shape as mhs
    padded_mdp = np.zeros(larger_shape)  # Assuming mdp has the same shape as mhs
    padded_mwn = np.zeros(larger_shape)  # Assuming mwn has the same shape as mhs

    padded_btm = np.zeros(larger_shape)  # Assuming btm has the same shape as bhs
    padded_btp = np.zeros(larger_shape)  # Assuming btp has the same shape as bhs
    padded_bdm = np.zeros(larger_shape)  # Assuming bdm has the same shape as bhs
    padded_bdp = np.zeros(larger_shape)  # Assuming bdp has the same shape as bhs
    padded_bwind = np.zeros(larger_shape)  # Assuming bwind has the same shape as bhs

    # Copy values from original arrays to padded arrays for other variables
    padded_mtm[:mtm.shape[0], :mtm.shape[1]] = mtm
    padded_mtp[:mtp.shape[0], :mtp.shape[1]] = mtp
    padded_mdm[:mdm.shape[0], :mdm.shape[1]] = mdm
    padded_mdp[:mdp.shape[0], :mdp.shape[1]] = mdp
    padded_mwn[:mwn.shape[0], :mwn.shape[1]] = mwn

    padded_btm[:btm.shape[0], :btm.shape[1]] = btm
    padded_btp[:btp.shape[0], :btp.shape[1]] = btp
    padded_bdm[:bdm.shape[0], :bdm.shape[1]] = bdm
    padded_bdp[:bdp.shape[0], :bdp.shape[1]] = bdp
    padded_bwind[:bwind.shape[0], :bwind.shape[1]] = bwind



#	ind=np.where( (mhs>0.0) & (bhs>0.0) )
    ind = np.where((padded_mhs > 0.0) & (padded_bhs > 0.0))
if np.size(ind)>0:
	print(' Total amount of matchups model/buoy: '+repr(np.size(ind)))

	# Save netcdf output file 
	lon[lon>180.]=lon[lon>180.]-360.
	initime=str(time.gmtime(mtime.min())[0])+str(time.gmtime(mtime.min())[1]).zfill(2)+str(time.gmtime(mtime.min())[2]).zfill(2)+str(time.gmtime(mtime.min())[3]).zfill(2)
	fintime=str(time.gmtime(mtime.max())[0])+str(time.gmtime(mtime.max())[1]).zfill(2)+str(time.gmtime(mtime.max())[2]).zfill(2)+str(time.gmtime(mtime.max())[3]).zfill(2)
	ncfile = nc.Dataset('WW3.Buoy'+ftag+'_'+initime+'to'+fintime+'.nc', "w", format=fnetcdf)
	ncfile.history="Matchups of WAVEWATCHIII point output (table) and NDBC and Copernicus Buoys. Total of "+repr(bhs[bhs>0.].shape[0])+" observations or pairs model/observation."
	# create  dimensions
	ncfile.createDimension('buoypoints', bhs.shape[0] )
	if gridinfo!=0:
		ncfile.createDimension('GlobalOceansSeas', ocnames.shape[0] )
		ncfile.createDimension('HighSeasMarineZones', hsmznames.shape[0] )
	if cyclonemap!=0:
		ncfile.createDimension('cycloneinfo', cinfo.shape[0] )
		vcinfo = ncfile.createVariable('cycloneinfo',dtype('a25'),('cycloneinfo'))
	# create variables.
	vstname = ncfile.createVariable('buoyID',dtype('a25'),('buoypoints'))
	vlat = ncfile.createVariable('latitude',np.dtype('float32').char,('buoypoints'))
	vlon = ncfile.createVariable('longitude',np.dtype('float32').char,('buoypoints'))

	if forecastds>0:
		ncfile.createDimension('time', nmhs.shape[2] )
		ncfile.createDimension('fcycle', unt.shape[0] )
		vt = ncfile.createVariable('time',np.dtype('float64').char,('fcycle','time'))
		vmhs = ncfile.createVariable('model_hs',np.dtype('float32').char,('buoypoints','fcycle','time'))
		vmtm = ncfile.createVariable('model_tm',np.dtype('float32').char,('buoypoints','fcycle','time'))
		vmtp = ncfile.createVariable('model_tp',np.dtype('float32').char,('buoypoints','fcycle','time'))
		vmdm = ncfile.createVariable('model_dm',np.dtype('float32').char,('buoypoints','fcycle','time'))
		vmdp = ncfile.createVariable('model_dp',np.dtype('float32').char,('buoypoints','fcycle','time'))
		vmwn = ncfile.createVariable('model_wind',np.dtype('float32').char,('buoypoints','fcycle','time'))

		vbhs = ncfile.createVariable('obs_hs',np.dtype('float32').char,('buoypoints','fcycle','time'))
		vbtm = ncfile.createVariable('obs_tm',np.dtype('float32').char,('buoypoints','fcycle','time'))
		vbtp = ncfile.createVariable('obs_tp',np.dtype('float32').char,('buoypoints','fcycle','time'))
		vbdm = ncfile.createVariable('obs_dm',np.dtype('float32').char,('buoypoints','fcycle','time'))
		vbdp = ncfile.createVariable('obs_dp',np.dtype('float32').char,('buoypoints','fcycle','time'))
		vbwind = ncfile.createVariable('obs_wind',np.dtype('float32').char,('buoypoints','fcycle','time'))

	else:
		ncfile.createDimension('time', bhs.shape[1] )
		vt = ncfile.createVariable('time',np.dtype('float64').char,('time'))
		vmhs = ncfile.createVariable('model_hs',np.dtype('float32').char,('buoypoints','time'))
		vmtm = ncfile.createVariable('model_tm',np.dtype('float32').char,('buoypoints','time'))
		vmtp = ncfile.createVariable('model_tp',np.dtype('float32').char,('buoypoints','time'))
		vmdm = ncfile.createVariable('model_dm',np.dtype('float32').char,('buoypoints','time'))
		vmdp = ncfile.createVariable('model_dp',np.dtype('float32').char,('buoypoints','time'))
		vmwn = ncfile.createVariable('model_wind',np.dtype('float32').char,('buoypoints','time'))

		vbhs = ncfile.createVariable('obs_hs',np.dtype('float32').char,('buoypoints','time'))
		vbtm = ncfile.createVariable('obs_tm',np.dtype('float32').char,('buoypoints','time'))
		vbtp = ncfile.createVariable('obs_tp',np.dtype('float32').char,('buoypoints','time'))
		vbdm = ncfile.createVariable('obs_dm',np.dtype('float32').char,('buoypoints','time'))
		vbdp = ncfile.createVariable('obs_dp',np.dtype('float32').char,('buoypoints','time'))
		vbwind = ncfile.createVariable('obs_wind',np.dtype('float32').char,('buoypoints','time'))




	if gridinfo!=0:
		vpdistcoast = ncfile.createVariable('distcoast',np.dtype('float32').char,('buoypoints')) 
		vpdepth = ncfile.createVariable('depth',np.dtype('float32').char,('buoypoints')) 
		vponi = ncfile.createVariable('GlobalOceansSeas',np.dtype('float32').char,('buoypoints'))
		vocnames = ncfile.createVariable('names_GlobalOceansSeas',dtype('a25'),('GlobalOceansSeas'))
		vphsmz = ncfile.createVariable('HighSeasMarineZones',np.dtype('float32').char,('buoypoints')) 
		vhsmznames = ncfile.createVariable('names_HighSeasMarineZones',dtype('a25'),('HighSeasMarineZones'))
	if cyclonemap!=0:
		if forecastds>0:
			vcmap = ncfile.createVariable('cyclone',np.dtype('float32').char,('buoypoints','fcycle','time'))
		else:
			vcmap = ncfile.createVariable('cyclone',np.dtype('float32').char,('buoypoints','time'))

	# Assign units
	vlat.units = 'degrees_north' ; vlon.units = 'degrees_east'
	vt.units = 'seconds since 1970-01-01T00:00:00+00:00'
	vmhs.units='m'; vbhs.units='m'
	vmtm.units='s'; vbtm.units='s'
	vmtp.units='s'; vbtp.units='s'
	vmdm.units='degrees'; vbdm.units='degrees'
	vmdp.units='degrees'; vbdp.units='degrees'
	vmwn.unit='m/s';vbwind.unit='m/s'

	if gridinfo!=0:
		vpdepth.units='m'; vpdistcoast.units='km'

	# Allocate Data
	vstname[:]=stname[:]; vlat[:] = lat[:]; vlon[:] = lon[:]
	if forecastds>0:
		vt[:,:]=nmtime[:,:]
		vmhs[:,:,:]=nmhs[:,:,:]
		vmtm[:,:,:]=nmtm[:,:,:]
		vmtp[:,:,:]=nmtp[:,:,:]
		vmdm[:,:,:]=nmdm[:,:,:]
		vmdp[:,:,:]=nmdp[:,:,:]
		vmwn[:,:,:]=nmwn[:,:,:]
		vbhs[:,:,:]=nbhs[:,:,:]
		vbtm[:,:,:]=nbtm[:,:,:]
		vbtp[:,:,:]=nbtp[:,:,:]
		vbdm[:,:,:]=nbdm[:,:,:]
		vbdp[:,:,:]=nbdp[:,:,:]
		vbwind[:,:,:]=nbwind[:,:,:]

	else:
		vt[:]=mtime[:]
		vmhs[:,:]=padded_mhs[:,:]
		vmtm[:,:]=padded_mtm[:,:]
		vmtp[:,:]=padded_mtp[:,:]
		vmdm[:,:]=padded_mdm[:,:]
		vmdp[:,:]=padded_mdp[:,:]
		vmwn[:,:]=padded_mwn[:,:]
		vbhs[:,:]=padded_bhs[:,:]
		vbtm[:,:]=btm[:,:]
		vbtp[:,:]=btp[:,:]
		vbdm[:,:]=bdm[:,:]
		vbdp[:,:]=bdp[:,:]
		vbwind[:,:]=bwind[:,:]




	if gridinfo!=0:
		vpdistcoast[:]=pdistcoast[:]
		vpdepth[:]=pdepth[:]
		vponi[:]=poni[:]; vocnames[:] = ocnames[:]
		vphsmz[:]=phsmz[:]; vhsmznames[:] = hsmznames[:]
	if cyclonemap!=0:
		vcinfo[:] = cinfo[:]
		if forecastds>0:
			vcmap[:,:,:]=nfcmap[:,:,:]	
		else:
			vcmap[:,:]=fcmap[:,:]

	ncfile.close()
	print(' ')
	print('Done. Netcdf ok. New file saved: WW3.Buoy'+ftag+'_'+initime+'to'+fintime+'.nc')
