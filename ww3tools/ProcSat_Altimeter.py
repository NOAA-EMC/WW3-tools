import numpy as np
from matplotlib.mlab import *
from pylab import *
import pandas as pd
import os
import netCDF4 as nc
import xarray as xr
import pyresample
from mpl_toolkits.basemap import shiftgrid
from geopy.distance import geodesic
import time
import timeit
from datetime import datetime
from calendar import timegm
import wread
import sys
import warnings; warnings.filterwarnings("ignore")


def along_track(AODN,atime,wconfig):
    '''
    The target is the observation space.
    Processing averages along-track, holding the obs space as reference.
    Inputs:
     (1) AODN pandas dataframe from wread.aodn_altimeter
     (2) time array (seconds since 1970)
     (3) wconfig dictionary, from wread.readconfig('ww3tools.yaml')
    Output: pandas dataframe containing: TIME (seconds since 1970), SAT_TIME (seconds since 1970),
      LATITUDE, LONGITUDE, WDEPTH, DISTCOAST, HS, HS_CAL, WSPD, WSPD_CAL
    Maryam Mohammadpour & Ricardo M. Campos
    '''

    # JASON3  112 minutes per revolution around the earth
    # https://www.ecmwf.int/en/newsletter/149/meteorology/use-radar-altimeter-products-ecmwf
    # https://www.nesdis.noaa.gov/current-satellite-missions/currently-flying/jason-3/jason-3-mission







    # spatial averages along track, using wconfig['dlim']
    i=0
    stime=[]; rstime=[]; slat=[]; slon=[]; swdepth=[]; swdistcoast=[]
    shs=[]; shscal=[]; swsp=[]; swspcal=[]
    for t in range(0,len(atime)):


        indt = np.where( abs(AODN['TIME'][:]-atime[t]) < wconfig['maxti'] )


        if np.size(indt)>1:
            indt=indt[0]
            # plot(AODN['LATITUDE'][indt],AODN['LONGITUDE'][indt],'.')
            lat=np.array(AODN['LATITUDE'][indt])
            lon=np.array(AODN['LONGITUDE'][indt]); lon[lon>360]=lon[lon>360]-360.

            # Coordinates of a point
            firstcoord = (lat[0], lon[0])
            segtrack = list(zip(lat, lon))
            # Calculate distances to reference coordinates
            distances = np.array([geodesic(firstcoord, ref_coords).meters for ref_coords in segtrack])
            # plot(distances,'.')
            segments = np.array(distances // wconfig['dlim']).astype('int')
            # print(' Total of '+repr(len(np.unique(segments)))+' segments')
            # 1h (+-30min) JASON3 2.850 Kms (487 records). Total of 82 segments of 35km
            for j in range(0,len(np.unique(segments))):
                indseg=np.where(segments==int(j))
                if np.size(indseg)>=1:
                    indseg=indseg[0]

                    stime = np.append(stime,atime[t])
                    rstime = np.append(rstime,np.nanmean(AODN['TIME'][indt[indseg]]))
                    slat = np.append(slat,np.nanmean(AODN['LATITUDE'][indt[indseg]]))
                    slon = np.append(slon,np.array(AODN['LONGITUDE'][indt[indseg]])[int(np.round(len(indseg)/2))]) # avoid problems with 360-0
                    swdepth = np.append(swdepth,np.nanmean(AODN['WDEPTH'][indt[indseg]]))
                    swdistcoast = np.append(swdistcoast,np.nanmean(AODN['DISTCOAST'][indt[indseg]]))
                    shs = np.append(shs,np.nanmean(AODN['HS'][indt[indseg]]))
                    shscal = np.append(shscal,np.nanmean(AODN['HS_CAL'][indt[indseg]]))
                    swsp = np.append(swsp,np.nanmean(AODN['WSPD'][indt[indseg]]))
                    swspcal = np.append(swspcal,np.nanmean(AODN['WSPD_CAL'][indt[indseg]]))

                    # print(repr(t)+" "+repr(i))
                    i=i+1
    swdepth = np.array(swdepth)
    swdistcoast = np.array(swdistcoast)
    shs = np.array(shs)
    swsp = np.array(swsp)

    if wconfig['qc']==0:
        indq = np.where( (stime>=adatemin) & (stime<=adatemax) )
    else:
        indq = np.where( (swdepth>=wconfig['mindepth']) & (swdistcoast>=wconfig['mindfc']) &
            (shs>wconfig['hsmin']) & (shs<wconfig['hsmax']) & (swsp>wconfig['wspmin']) & (swsp<wconfig['wspmax']) &
            (stime>=adatemin) & (stime<=adatemax) )

    if np.size(indq)>2:
        indq=indq[0]
        procaodn = {'TIME': stime[indq],'SAT_TIME': rstime[indq], 'LATITUDE': slat[indq], 'LONGITUDE': slon[indq],
            'WDEPTH': swdepth[indq], 'DISTCOAST': swdistcoast[indq],
            'HS': shs[indq], 'HS_CAL': shscal[indq],
            'WSPD': swsp[indq], 'WSPD_CAL': swspcal[indq]}

    else:
        procaodn = {'TIME': [],'SAT_TIME': [], 'LATITUDE': [], 'LONGITUDE': [],
            'WDEPTH': [], 'DISTCOAST': [],
            'HS': [], 'HS_CAL': [],
            'WSPD': [], 'WSPD_CAL': []}

    AODN_ALONGTRACK = pd.DataFrame(procaodn)
    return AODN_ALONGTRACK

# weight function for pyresample
def wf(pdist):
    a=(1 - pdist / (wconfig['dlim']+1))
    return (abs(a)+a)/2

def gridded(AODN,atime,wconfig):
    '''
    The target is the model space.
    Processing satellite collocation into grid points (structured or unstructured).
    Inputs:
     (1) AODN pandas dataframe from wread.aodn_altimeter
     (2) time array (seconds since 1970)
     (3) wconfig dictionary, from wread.readconfig('ww3tools.yaml')
    Output: pandas dataframe containing: TIME (seconds since 1970), SAT_TIME (seconds since 1970),
      LATITUDE, LONGITUDE, WDEPTH, DISTCOAST, HS, HS_CAL, WSPD, WSPD_CAL
    Ricardo M. Campos
    '''

    if str(wconfig['grid_info']).split('.')[-1] == 'grib2' or str(wconfig['grid_info']).split('.')[-1] == 'grb2':
        # grib2 format
        ds = xr.open_dataset(wconfig['grid_info'], engine='cfgrib')
    else:
        # netcdf format
        ds = xr.open_dataset(wconfig['grid_info'])

    if ('latitude' in list(ds.coords)) or ('latitude' in list(ds.keys())):
        lat = np.array(ds['latitude'].values)
    elif ('LATITUDE' in list(ds.coords)) or ('LATITUDE' in list(ds.keys())):
        lat = np.array(ds['LATITUDE'].values)
    elif ('lat' in list(ds.coords)) or ('lat' in list(ds.keys())):
        lat = np.array(ds['lat'].values)
    elif ('LAT' in list(ds.coords)) or ('LAT' in list(ds.keys())):
        lat = np.array(ds['LAT'].values)

    if ('longitude' in list(ds.coords)) or ('longitude' in list(ds.keys())):
        lon = np.array(ds['longitude'].values)
    elif ('LONGITUDE' in list(ds.coords)) or ('LONGITUDE' in list(ds.keys())):
        lon = np.array(ds['LONGITUDE'].values)
    elif ('lon' in list(ds.coords)) or ('lon' in list(ds.keys())):
        lon = np.array(ds['lon'].values)
    elif ('LON' in list(ds.coords)) or ('LON' in list(ds.keys())):
        lon = np.array(ds['LON'].values)

    lon[lon>180]=lon[lon>180]-360.
    AODN['LONGITUDE'][AODN['LONGITUDE']>180.]=AODN['LONGITUDE'][AODN['LONGITUDE']>180.]-360.

    if wconfig['grid_format']==1:
        # Structured grid
        NEWLON, NEWLAT = meshgrid(lon,lat)
        NEWLON = np.reshape(NEWLON,(len(lon)*len(lat)))
        NEWLAT = np.reshape(NEWLAT,(len(lon)*len(lat)))
    else:
        # Unstructured grid
        NEWLON, NEWLAT = lon,lat

    if 'mask' in list(ds.keys()):
        mask=np.array(ds['mask'].values)
        if wconfig['grid_format']==1:
            mask = np.reshape(mask,(mask.shape[0]*mask.shape[1]))

        ind=np.where(mask>0)[0]; flatm = np.array(NEWLAT[ind]); flonm = np.array(NEWLON[ind])
        del ind
    else:
        flatm = np.array(NEWLAT); flonm = np.array(NEWLON)

    del NEWLON, NEWLAT
    targ_def = pyresample.geometry.SwathDefinition(lons=flonm,lats=flatm)

    slen=len(AODN); fctr=int(np.ceil(wconfig['dlim']/(np.nanmean(np.diff(np.unique(flatm)))*(10**8)))**2)
    ftime=np.double(np.zeros((slen*fctr),'d')); flat=np.zeros((slen*fctr),'f'); flon=np.zeros((slen*fctr),'f')
    fhsk=np.zeros((slen*fctr),'f'); stdhsk=np.zeros((slen*fctr),'f'); counthsk=np.zeros((slen*fctr),'f')
    fhskcal=np.zeros((slen*fctr),'f'); fwnd=np.zeros((slen*fctr),'f'); fwndcal=np.zeros((slen*fctr),'f')
    fwdepth=np.zeros((slen*fctr),'f'); fdfc=np.zeros((slen*fctr),'f');

    # into the regular grid with pyresample kd tree
    ii=0
    for t in range(0,len(atime)):
        indt = np.where( abs(AODN['TIME'].values[:]-atime[t]) < wconfig['maxti'] )
        if np.size(indt)>1:
            indt=indt[0]
            prlon=np.copy(AODN['LONGITUDE'].values[indt])
            orig_def = pyresample.geometry.SwathDefinition(lons=prlon, lats=AODN['LONGITUDE'].values[indt]); del prlon
            # By distance function wf
            auxfhsk, auxstdhsk, auxcounthsk = pyresample.kd_tree.resample_custom(orig_def,AODN['HS'].values[indt],targ_def,radius_of_influence=wconfig['dlim'],weight_funcs=wf,fill_value=0,with_uncert=True,nprocs=wconfig['npcs'])
            auxfhskcal = pyresample.kd_tree.resample_custom(orig_def,AODN['HS_CAL'].values[indt],targ_def,radius_of_influence=wconfig['dlim'],weight_funcs=wf,fill_value=0,nprocs=wconfig['npcs'])
            auxfwnd = pyresample.kd_tree.resample_custom(orig_def,AODN['WSPD'].values[indt],targ_def,radius_of_influence=wconfig['dlim'],weight_funcs=wf,fill_value=0,nprocs=wconfig['npcs'])
            auxfwndcal = pyresample.kd_tree.resample_custom(orig_def,AODN['WSPD_CAL'].values[indt],targ_def,radius_of_influence=wconfig['dlim'],weight_funcs=wf,fill_value=0,nprocs=wconfig['npcs'])
            auxfwdepth = pyresample.kd_tree.resample_nearest(orig_def,AODN['WDEPTH'].values[indt],targ_def,radius_of_influence=wconfig['dlim'],fill_value=0,nprocs=wconfig['npcs'])
            auxfdfc = pyresample.kd_tree.resample_nearest(orig_def,AODN['DISTCOAST'].values[indt],targ_def,radius_of_influence=wconfig['dlim'],fill_value=0,nprocs=wconfig['npcs'])
            # print('   - ok '+repr(t))
            indpqq = np.where( (auxfhskcal>0.01) & (auxfwnd>0.01) & (auxfhskcal<wconfig['hsmax']) & (auxfwnd<wconfig['wspmax']) &
                (auxstdhsk>0.001) & (auxcounthsk>0.) )

            if np.size(indpqq)>0:
                indpqq=indpqq[0]
                # allocate data into final array
                ftime[ii:ii+np.size(indpqq)] = np.array(np.zeros(np.size(indpqq),'d')+atime[t]).astype('double')
                flat[ii:ii+np.size(indpqq)] = np.array(flatm[indpqq]).astype('float')
                flon[ii:ii+np.size(indpqq)] = np.array(flonm[indpqq]).astype('float')
                fhsk[ii:ii+np.size(indpqq)] = np.array(auxfhsk[indpqq]).astype('float')
                fhskcal[ii:ii+np.size(indpqq)] = np.array(auxfhskcal[indpqq]).astype('float')
                fwnd[ii:ii+np.size(indpqq)] = np.array(auxfwnd[indpqq]).astype('float')
                fwndcal[ii:ii+np.size(indpqq)] = np.array(auxfwndcal[indpqq]).astype('float')
                fwdepth[ii:ii+np.size(indpqq)] = np.array(auxfwdepth[indpqq]).astype('float')
                fdfc[ii:ii+np.size(indpqq)] = np.array(auxfdfc[indpqq]).astype('float')
                ii=ii+len(indpqq)

            del auxfhsk, auxstdhsk, auxcounthsk, auxfhskcal, auxfwnd, auxfwndcal, indpqq

        del indt
        print('PyResample kdtree, hourly time, '+repr(t)+' of '+repr(atime.shape[0]))

    if wconfig['qc']==0:
        indf=np.where( (ftime>0.) & (fhsk>=0.0) )
    else:
        ind=np.where( (fhsk<wconfig['hsmin']) | (fhsk>wconfig['hsmax']) )
        if np.any(ind):
            fhsk[ind[0]]=np.nan; del ind

        ind=np.where( (fhskcal<wconfig['hsmin']) | (fhskcal>wconfig['hsmax']) )
        if np.any(ind):
            fhskcal[ind[0]]=np.nan; del ind

        ind=np.where( (fwnd<wconfig['wspmin']) | (fwnd>wconfig['wspmax']) )
        if np.any(ind):
            fwnd[ind[0]]=np.nan; del ind

        ind=np.where( (fwndcal<wconfig['wspmin']) | (fwndcal>wconfig['wspmax']) )
        if np.any(ind):
            fwndcal[ind[0]]=np.nan; del ind

        indf=np.where( (ftime>0.) & (fhsk>=wconfig['hsmin']) & (fwdepth>=wconfig['mindepth']) & (fdfc>=wconfig['mindfc']) )

    if np.size(indf)>0:
        indf=indf[0]
        ftime=np.array(ftime[indf]).astype('double')
        flat=np.array(flat[indf]); flon=np.array(flon[indf])
        fwdepth=np.array(fwdepth[indf]); fdfc=np.array(fdfc[indf])
        fhsk=np.array(fhsk[indf]); fhskcal=np.array(fhskcal[indf])
        fwnd=np.array(fwnd[indf]); fwndcal=np.array(fwndcal[indf])

        procaodn = {'TIME': ftime,'SAT_TIME': ftime, 'LATITUDE': flat, 'LONGITUDE': flon,
            'WDEPTH': fwdepth, 'DISTCOAST': fdfc,
            'HS': fhsk, 'HS_CAL': fhskcal,
            'WSPD': fwnd, 'WSPD_CAL': fwndcal}

    else:
        procaodn = {'TIME': [],'SAT_TIME': [], 'LATITUDE': [], 'LONGITUDE': [],
            'WDEPTH': [], 'DISTCOAST': [],
            'HS': [], 'HS_CAL': [],
            'WSPD': [], 'WSPD_CAL': []}

    print(' '); print(' gridded() Competed. Array Size '+str(np.size(indf))); print(' ')
    AODN_GRIDDED = pd.DataFrame(procaodn)
    return AODN_GRIDDED


def savesat(AODN,wconfig,altsel):
    '''
    Save AODN processed data (along-track or gridded)
    Inputs:
     (1) AODN pandas dataframe
     (2) wconfig dictionary, from wread.readconfig('ww3tools.yaml')
     (3) altimeter name
    Output: netcdf and csv (if array is not huge) saved in the output path
     defined in ww3tools.yaml
    '''

    if wconfig['tspace']==1:
        smethod="Gridded"
        hmsg="AODN Altimeter data on model grid."
    elif wconfig['tspace']==2:
        smethod="AlongTrack"
        hmsg="AODN Altimeter data processed along track."

    datein = datetime.utcfromtimestamp(AODN_ALONGTRACK['TIME'][0]).strftime('%Y%m%d%H')
    datefin = datetime.utcfromtimestamp(AODN_ALONGTRACK['TIME'].iloc[-1]).strftime('%Y%m%d%H')

    fname=wconfig['path_out']+"Altimeter"+smethod+"_"+wconfig['ftag']+"_"+altsel+"_"+datein+"to"+datefin

    # Save netcdf
    ncfile = nc.Dataset(fname+".nc", "w", format=wconfig['fnetcdf'])
    ncfile.history=hmsg
    # create  dimensions.
    ncfile.createDimension('time' , len(AODN['TIME']))
    ncfile.createDimension('sname' , 1 )
    # create variables.
    vflat = ncfile.createVariable('latitude',np.dtype('float32').char,('time'))
    vflon = ncfile.createVariable('longitude',np.dtype('float32').char,('time'))
    vft = ncfile.createVariable('time',np.dtype('float64').char,('time'))
    vfst = ncfile.createVariable('sat_time',np.dtype('float64').char,('time'))
    vsname = ncfile.createVariable('sat_name',dtype('a25'),('sname'))
    vfhs = ncfile.createVariable('hs',np.dtype('float32').char,('time'))
    vfhscal = ncfile.createVariable('hs_cal',np.dtype('float32').char,('time'))
    vfwnd = ncfile.createVariable('wsp',np.dtype('float32').char,('time'))
    vfwndcal = ncfile.createVariable('wsp_cal',np.dtype('float32').char,('time'))
    # Assign units
    vft.units = 'seconds since 1970-01-01 00:00:00'
    vfst.units = 'seconds since 1970-01-01 00:00:00'
    vflat.units = 'degrees_north' ; vflon.units = 'degrees_east'
    vfwnd.units = 'm/s' ; vfwndcal.units = 'm/s'
    vfhs.units = 'm'; vfhscal.units = 'm'
    # Allocate Data
    vsname[:] = np.array(altsel).astype('str')
    vflat[:] = np.array(AODN['LATITUDE']); vflon[:] = np.array(AODN['LONGITUDE'])
    vft[:] = np.array(AODN['TIME']).astype('double'); vfst[:] = np.array(AODN['SAT_TIME']).astype('double');
    vfhs[:] = np.array(AODN['HS']); vfhscal[:] = np.array(AODN['HS_CAL'])
    vfwnd[:] = np.array(AODN['WSPD']); vfwndcal[:] = np.array(AODN['WSPD_CAL'])
    ncfile.close()
    print(' ')
    print("netcdf ok: "+fname+".nc")

    del vft, vfst, vflat, vflon, vfhs,vfhscal, vfwnd, vfwndcal

    if len(AODN['TIME'])<10**8 :
        AODN.iloc[:, 2:] = AODN.iloc[:, 2:].round(3)
        AODN.to_csv(fname+".csv", index=False)

    print("text file ok: "+fname+".csv")


if __name__ == "__main__":

    # start time
    start = timeit.default_timer()
    # WW3-tools configuration file
    wconfig=wread.readconfig('ww3tools.yaml')

    # Default date interval
    datemin='1985010100'; datemax = datetime.utcnow().strftime("%Y%m%d%H")
    # Default time step (in seconds) to build the final time array (regular)
    time_step=float(3600.)
    ww3grid=None
    # Read input arguments
    if len(sys.argv) < 2 :
        sys.exit(' At least 1 argument (satellite mission) must be entered.')
    if len(sys.argv) >= 2 :
        altsel=str(sys.argv[1]) # selection of the altimeter mission, see sdname below
    if len(sys.argv) >= 3 :
        datemin=str(sys.argv[2]) # initial date (YYYYMMDDHH)
    if len(sys.argv) >= 4 :
        datemax=str(sys.argv[3]) # final date (YYYYMMDDHH)
    if len(sys.argv) >= 5 :
        time_step=float(sys.argv[4]) # time step (in seconds) to build the final time array (regular)

    # read and organize AODN altimeter data for the mission and domain of interest
    AODN = wread.aodn_altimeter(altsel,wconfig,datemin,datemax)
    # AODN.to_csv(wconfig['path_out']+"AODN_altimeterSelection_"+datemin+"to"+datemax+".csv", index=False)

    # date interval in seconds since 1970, user selection
    adatemin= np.double(timegm( time.strptime(datemin, '%Y%m%d%H')))
    adatemax= np.double(timegm( time.strptime(datemax, '%Y%m%d%H')))
    atime = np.array(np.arange(adatemin,adatemax+1.,time_step)).astype('double')

    if wconfig['tspace']==2:
        AODN_ALONGTRACK = along_track(AODN,atime,wconfig)
        savesat(AODN_ALONGTRACK,wconfig,altsel)

    elif wconfig['tspace']==1:
        AODN_GRID = gridded(AODN,atime,wconfig)
        savesat(AODN_GRID,wconfig,altsel)

    stop = timeit.default_timer()
    print('ProcSat_Altimeter.py  successfully completed in '+repr(int(round(stop - start,0)))+' seconds. See the output at '+wconfig['path_out'])



