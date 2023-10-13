#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
wproc.py

VERSION AND LAST UPDATE:
 v1.0  07/10/2023
 v1.1  07/24/2023
 v1.2  10/11/2023

PURPOSE:
 Auxiliary processing functions.
 Functions:
   orgensemblesat

USAGE and OUTPUT:
 The explanation for each function is contained in the headers

DEPENDENCIES:
 See setup.py and the imports below.

AUTHOR and DATE:
 07/10/2023: Ricardo M. Campos, first version.
 07/24/2023: Ricardo M. Campos, debug orgensemblesat and improve processing speed. 
 10/11/2023: Ricardo M. Campos, new functions to support modularization.

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import os
import netCDF4 as nc
import numpy as np
import sys
import pandas as pd
import os
import warnings; warnings.filterwarnings("ignore")


def wlevconv(data_lev1=None,lev1=None,lev2=None,alpha=0.12):
    """
    Wind height conversion DNV-C209 standard/recommendation
    Input:
    - time-series (numpy array) of origin level
    - height (meters) of origin level
    - height (meters) of target level
    - alpha (see DNV-C209 standard/recommendation)
    Output:
    - time-series converted to the target level

    Example1 (4.1 m height obs data to 10-m model level):
    lev1=4.1; lev2=10; data_lev1=10 #(m/s)
    converted_wind = wlevconv(data_lev1,lev1,lev2)

    Example2 (10-m wind from model to 3.8 m anemometer position):
    lev1=10; lev2=3.8; data_lev1=10 #(m/s)
    converted_wind = wlevconv(data_lev1,lev1,lev2)
    """

    if np.any(lev1)==None or np.any(lev2)==None or np.any(data_lev1)==None:
        raise ValueError("Two levels (in meters) and one data record (wind in m/s) must be informed.")

    mfactor=((lev2/lev1)**(alpha))
    data_lev2 = (mfactor * data_lev1)

    return data_lev2
    print(' Wind conversion ok')


def interp_nan(data,lmt=10**50):
    '''
    Fill NaN values with linear interpolation.
    User can enter one or two inputs:
      1) time-series containing NaN values to fill in
      2) maximum number of consecutive NaN values to interpolate (to
         avoid interpolating long segments)
    '''
    data=np.array(data)
    lmt=int(lmt)

    if data.ndim>1:
        raise ValueError(' Input array with too many dimensions. Only time-series (1 dimension) allowed.')
    else:
        indd=np.where(data>-999.)
        if np.size(indd)>0:
            indd=indd[0][-1] # last valid number
            adata=np.array(data[0:indd])

            # using pandas
            A=pd.Series(adata)
            # B=A.interpolate(method="polynomial",order=2,limit=lmt)
            B=A.interpolate(method="linear",limit=lmt)
            B=np.array(B.values)

            data[0:indd]=np.array(B)
           
            del lmt,A,B

    return data


def orgensemblesat(flist,nmb,esize='yes'):
    '''
    This function reads and organizes the ensemble members and 
    extra information generated from multiple runs of 
    modelSat_collocation.py (one for each ensemble member).
    The netcdf files are for example:
    WW3.Altimeter_200001_c00_2000011203to2000030100.nc
    where c00 to p10 are the ensemble members. You can also enter
    the full path (recommended).
    This function reads a list (string array or text file) as the 
    first argument containing the files (can be only one) for one 
    specific single member.
    The code will read the file string and loop through the ensemble
    members (the number of members is the second input argument).
    The format must be WW3.Altimeter_DATE_c00 for the control member
    and WW3.Altimeter_DATE_p01 for all the perturbed members, saved
    in the same directory. The split string separates using "_" so
    any date format can be entered.
    You can create a list with:
    ls -d $PWD/*_p10_*.nc > list.txt &
    and input the list.txt as the first argument,
    or read 'list.txt' with flist = np.atleast_1d(np.loadtxt('list.txt',dtype=str))
    and enter flist as the first argument.
    Inputs:
     (1) list of WW3.Altimeter files for one specific member.
     (2) number of members. Remember to count it properly, including 0 as the control member.
    Output:
     (1) One dataframe containing the control member (ONLY) and extra
      information, such as lat, lon, distcoast, depth etc. The shape[1]
      of this dataframe contain the extra information and NOT the ensemble
      members.
     (2) Two numpy arrays (Hs and U10) with all members included into
      a new dimension (columns, shape[1]).
     (3) Dictionary with names of GlobalOceansSeas, HighSeasMarineZones,
      satellites, and cycloneinfo.
    Example:

     import wproc
     fdatacontrol,fmhs,fmwnd,dfnames=wproc.orgensemblesat('list.txt',11)

      you can see the information available with fdatacontrol.keys()
      where fdatacontrol['obs_hs'].values is the array with satellite observations of Hs
      and fdatacontrol['model_hs'].values is the array with model results of Hs
      for the control member.
     You can directly compose the whole ensemble (including control and
      perturbed members) with: np.c_[fdatacontrol['model_hs'].values,fmhs] where
      the shape[1] dimension contains the ensemble member.
      The ensemble mean will be np.nanmean(np.c_[fdatacontrol['model_hs'].values,fmhs],axis=1)

     If you want to save fdatacontrol,fmhs,fmwnd for latter use, you can do using pickle format:
     import pickle
     with open("WW3Altimeter_Data.pickle", 'wb') as file:
      pickle.dump(fdatacontrol, file)
      pickle.dump(fmhs, file)
      pickle.dump(fmwnd, file)
     and read later on using:
     with open('WW3Altimeter_Data.pickle', 'rb') as file:
      fdatacontrol = pickle.load(file)
      fmhs = pickle.load(file)
      fmwnd = pickle.load(file)
    '''

    # Check how flist is entered. It can be 'list.txt' directly or a numpy array containing the list.
    if isinstance(flist, str):
        flist = np.atleast_1d(np.loadtxt(flist,dtype=str))
    elif isinstance(flist, np.ndarray):
        flist = np.array(flist).astype('str')

    print(" "); print(" Function wproc.orgensemblesat. Read and organize data of "+repr(nmb)+" ensemble members ..."); print(" ")

    # ensemble members
    ensm=np.arange(0,nmb).astype('int')

    c=0
    for i in range(0,np.size(flist)):

        cpnd=1

        # Main loop: Ensemble members with exact match and equal array size.
        for j in ensm:

            if cpnd==1:
                if '/' in flist[i]:
                    if j==0:
                        os.system("ln -fs "+"/".join(flist[i].split('/')[0:-1])+"/"+"_".join(flist[i].split('/')[-1].split('_')[0:2])+"_c"+str(j).zfill(2)+"_*.nc"+" auxorgensemblesat.nc")
                    else:
                        os.system("ln -fs "+"/".join(flist[i].split('/')[0:-1])+"/"+"_".join(flist[i].split('/')[-1].split('_')[0:2])+"_p"+str(j).zfill(2)+"_*.nc"+" auxorgensemblesat.nc")
                else:
                    if j==0:
                        os.system("ln -fs "+"_".join(flist[i].split('_')[0:2])+"_c"+str(j).zfill(2)+"_*.nc"+" auxorgensemblesat.nc")
                    else:
                        os.system("ln -fs "+"_".join(flist[i].split('_')[0:2])+"_p"+str(j).zfill(2)+"_*.nc"+" auxorgensemblesat.nc")

                f=nc.Dataset("auxorgensemblesat.nc"); os.system("rm -f auxorgensemblesat.nc")

                if j==0:
                    df = pd.DataFrame({'time': np.array(f.variables['time'][:]).astype('double'),
                        'cycle': np.array(f.variables['cycle'][:]).astype('double'),
                        'latitude': np.array(f.variables['latitude'][:]),
                        'longitude': np.array(f.variables['longitude'][:]),
                        'satelliteID': np.array(f.variables['satelliteID'][:]),
                        'obs_hs': np.array(f.variables['obs_hs'][:]),
                        'obs_wnd': np.array(f.variables['obs_wnd'][:]),
                        'model_hs': np.array(f.variables['model_hs'][:]),
                        'model_wnd': np.array(f.variables['model_wnd'][:]),
                        'month': np.array(f.variables['month'][:]),
                        'distcoast': np.array(f.variables['distcoast'][:]),
                        'depth': np.array(f.variables['depth'][:]),
                        'GlobalOceansSeas': np.array(f.variables['GlobalOceansSeas'][:]),
                        'HighSeasMarineZones': np.array(f.variables['HighSeasMarineZones'][:]),
                        'cyclone': np.array(f.variables['cyclone'][:])},
                        columns=['time','cycle','latitude','longitude','satelliteID','obs_hs','obs_wnd','model_hs','model_wnd','month','distcoast','depth','GlobalOceansSeas','HighSeasMarineZones','cyclone'])

                    if c==0:
                        ocnames=np.array(f.variables['names_GlobalOceansSeas'][:]).astype('str')
                        hsmznames=np.array(f.variables['names_HighSeasMarineZones'][:]).astype('str')
                        sdname=np.array(f.variables['names_satellite'][:]).astype('str')
                        cinfo=np.array(f.variables['cycloneinfo'][:]).astype('str')
                        dfnames={'names_GlobalOceansSeas':list(ocnames),'names_HighSeasMarineZones':list(hsmznames),'names_satellite':list(sdname),'cycloneinfo':list(cinfo)}

                elif j==1:
                    if df.shape[0]==f.variables['model_hs'].shape[0]:
                        amhs=np.atleast_2d(np.array(f.variables['model_hs'][:]))
                        amwnd=np.atleast_2d(np.array(f.variables['model_wnd'][:]))
                    else:
                        print("  The size of arrays related to ensemble members (member "+repr(j)+") is different: "+"_".join(flist[i].split('_')[0:2]).split('/')[-1]+". Using pandas dataframe and pd.merge to select and organize the coincident records among different members.")
                        cpnd=0

                else:
                    if amhs.shape[1]==f.variables['model_hs'].shape[0]:
                        amhs=np.append(amhs,np.atleast_2d(np.array(f.variables['model_hs'][:])),axis=0)
                        amwnd=np.append(amwnd,np.atleast_2d(np.array(f.variables['model_wnd'][:])),axis=0)
                    else:
                        print("  The size of arrays related to ensemble members (member "+repr(j)+") is different: "+"_".join(flist[i].split('_')[0:2]).split('/')[-1]+". Using pandas dataframe and pd.merge to select and organize the coincident records among different members.")
                        cpnd=0
                        del amhs,amwnd,df

                f.close(); del f
                c=c+1

        # -------------------------------------------

        if cpnd==0:
            # Alternative loop: Ensemble members with different array sizes. -------------------------------------------
            print("  Alternative loop (ensemble members with different array sizes) "+"_".join(flist[i].split('_')[0:2]).split('/')[-1]+" ... ")
            # First build arrays of lat,lon,time,forecast time, and cycle time, common to every member:
            k=0
            for j in np.append(ensm,ensm[0]):
                if '/' in flist[i]:
                    if j==0:
                        os.system("ln -fs "+"/".join(flist[i].split('/')[0:-1])+"/"+"_".join(flist[i].split('/')[-1].split('_')[0:2])+"_c"+str(j).zfill(2)+"_*.nc"+" auxorgensemblesat.nc")
                    else:
                        os.system("ln -fs "+"/".join(flist[i].split('/')[0:-1])+"/"+"_".join(flist[i].split('/')[-1].split('_')[0:2])+"_p"+str(j).zfill(2)+"_*.nc"+" auxorgensemblesat.nc")
                else:
                    if j==0:
                        os.system("ln -fs "+"_".join(flist[i].split('_')[0:2])+"_c"+str(j).zfill(2)+"_*.nc"+" auxorgensemblesat.nc")
                    else:
                        os.system("ln -fs "+"_".join(flist[i].split('_')[0:2])+"_p"+str(j).zfill(2)+"_*.nc"+" auxorgensemblesat.nc")

                f=nc.Dataset("auxorgensemblesat.nc"); os.system("rm -f auxorgensemblesat.nc")
                if k==0:
                    dfa = pd.DataFrame({'time': np.array(f.variables['time'][:]).astype('double'),
                        'cycle': np.array(f.variables['cycle'][:]).astype('double'),
                        'latitude': np.array(f.variables['latitude'][:]).astype('double'),
                        'longitude': np.array(f.variables['longitude'][:]).astype('double'),
                        'satelliteID': np.array(f.variables['satelliteID'][:]).astype('double'),                    
                        'obs_hs': np.array(f.variables['obs_hs'][:]).astype('double')},
                        columns=['time','cycle','latitude','longitude','satelliteID','obs_hs'])

                    dfa = dfa[~dfa.duplicated(subset=['time','cycle','latitude','longitude','satelliteID','obs_hs'])]
                    dfa = dfa.reset_index(drop=True)

                dfb = pd.DataFrame({'time': np.array(f.variables['time'][:]).astype('double'),
                    'cycle': np.array(f.variables['cycle'][:]).astype('double'),
                    'latitude': np.array(f.variables['latitude'][:]).astype('double'),
                    'longitude': np.array(f.variables['longitude'][:]).astype('double'),
                    'satelliteID': np.array(f.variables['satelliteID'][:]).astype('double'),
                    'obs_hs': np.array(f.variables['obs_hs'][:]).astype('double')},
                    columns=['time','cycle','latitude','longitude','satelliteID','obs_hs'])

                dfb = dfb[~dfb.duplicated(subset=['time','cycle','latitude','longitude','satelliteID','obs_hs'])]
                dfb = dfb.reset_index(drop=True)

                dfa = pd.merge(dfa, dfb)
                dfa = dfa[~dfa.duplicated(subset=['time','cycle','latitude','longitude','satelliteID','obs_hs'])]
                dfa = dfa.reset_index(drop=True)
                f.close()
                del f,dfb; k=k+1

            # Now read and organize the different members using dfa as the reference of lat/lon/time/cycle/ftime
            for j in ensm:

                if '/' in flist[i]:
                    if j==0:
                        os.system("ln -fs "+"/".join(flist[i].split('/')[0:-1])+"/"+"_".join(flist[i].split('/')[-1].split('_')[0:2])+"_c"+str(j).zfill(2)+"_*.nc"+" auxorgensemblesat.nc")
                    else:
                        os.system("ln -fs "+"/".join(flist[i].split('/')[0:-1])+"/"+"_".join(flist[i].split('/')[-1].split('_')[0:2])+"_p"+str(j).zfill(2)+"_*.nc"+" auxorgensemblesat.nc")
                else:
                    if j==0:
                        os.system("ln -fs "+"_".join(flist[i].split('_')[0:2])+"_c"+str(j).zfill(2)+"_*.nc"+" auxorgensemblesat.nc")
                    else:
                        os.system("ln -fs "+"_".join(flist[i].split('_')[0:2])+"_p"+str(j).zfill(2)+"_*.nc"+" auxorgensemblesat.nc")

                f=nc.Dataset("auxorgensemblesat.nc"); os.system("rm -f auxorgensemblesat.nc")

                dfb = pd.DataFrame({'time': np.array(f.variables['time'][:]).astype('double'),
                    'cycle': np.array(f.variables['cycle'][:]).astype('double'),
                    'latitude': np.array(f.variables['latitude'][:]).astype('double'),
                    'longitude': np.array(f.variables['longitude'][:]).astype('double'),
                    'satelliteID': np.array(f.variables['satelliteID'][:]).astype('double'),
                    'obs_hs': np.array(f.variables['obs_hs'][:]).astype('double')},
                    columns=['time','cycle','latitude','longitude','satelliteID','obs_hs'])

                # fast way to find the elements of dfb contained in dfa
                dfb = pd.merge(dfb, dfa)
                ind = np.array(dfb.index).astype('int')
                del dfb

                if j==0:
                    df = pd.DataFrame({'time': np.array(f.variables['time'][:]).astype('double'),
                        'cycle': np.array(f.variables['cycle'][:]).astype('double'),
                        'latitude': np.array(f.variables['latitude'][:]),
                        'longitude': np.array(f.variables['longitude'][:]),
                        'satelliteID': np.array(f.variables['satelliteID'][:]),
                        'obs_hs': np.array(f.variables['obs_hs'][:]),
                        'obs_wnd': np.array(f.variables['obs_wnd'][:]),
                        'model_hs': np.array(f.variables['model_hs'][:]),
                        'model_wnd': np.array(f.variables['model_wnd'][:]),
                        'month': np.array(f.variables['month'][:]),
                        'distcoast': np.array(f.variables['distcoast'][:]),
                        'depth': np.array(f.variables['depth'][:]),
                        'GlobalOceansSeas': np.array(f.variables['GlobalOceansSeas'][:]),
                        'HighSeasMarineZones': np.array(f.variables['HighSeasMarineZones'][:]),
                        'cyclone': np.array(f.variables['cyclone'][:])},
                        columns=['time','cycle','latitude','longitude','satelliteID','obs_hs','obs_wnd','model_hs','model_wnd','month','distcoast','depth','GlobalOceansSeas','HighSeasMarineZones','cyclone'])

                    df = df.loc[ind]

                elif j==1:
                    amhs=np.atleast_2d(np.array(f.variables['model_hs'][:])[ind])
                    amwnd=np.atleast_2d(np.array(f.variables['model_wnd'][:])[ind])
                else:
                    amhs=np.append(amhs,np.atleast_2d(np.array(f.variables['model_hs'][:])[ind]),axis=0)
                    amwnd=np.append(amwnd,np.atleast_2d(np.array(f.variables['model_wnd'][:])[ind]),axis=0)

                f.close(); del ind,f

            del dfa
            # -------------------------------------------

        # fdata (control), and perturbed members: fmhs, fmwnd
        if i==0:
            fdata=df
            fmhs=np.copy(amhs.T)
            fmwnd=np.copy(amwnd.T)
        else:
            fdata = pd.concat([fdata, df]); fdata = fdata.reset_index(drop=True)
            fmhs=np.append(fmhs,amhs.T,axis=0)
            fmwnd=np.append(fmwnd,amwnd.T,axis=0)

        del df, amhs, amwnd
        print('  Ok  '+"_".join(flist[i].split('_')[0:2]).split('/')[-1]+'  for all ensemble members')

    print(" Done. wproc.orgensemblesat concluded.")
    return fdata,fmhs,fmwnd,dfnames
    # ------------------------


# Under development: spectral interpolation function, quality control



