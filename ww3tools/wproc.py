#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
wproc.py

VERSION AND LAST UPDATE:
 v1.0  07/10/2023

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

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import os
import netCDF4 as nc
import numpy as np
import sys
import pandas as pd
import warnings; warnings.filterwarnings("ignore")

def orgensemblesat(*args):
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
    in the save directory. The split string separates using "_" so
    any date format can be entered.
    You can create a list with:
    ls -d $PWD/*_p10_*.nc > list.txt &
    and input the list.txt as the first argument,
    or read the 'list.txt' with flist = np.atleast_1d(np.loadtxt('list.txt',dtype=str))
    and enter flist as the first argument.
    Inputs:
     (1) list of WW3.Altimeter files for one specific member.
     (2) number of members. Remember to count it properly, including 0 as the control member
    Output:
     (1) One dataframe containing the control member (ONLY) and extra
      information, such as lat, lon, distcoast, depth etc. The shape[1] 
      of this dataframe contain the extra information and NOT the ensemble
      members.
     (2) Two numpy arrays (Hs and U10) with all members included into
      a new dimension (columns, shape[1]).
    Example:
     import wproc
     fdatacontrol,fhs,fwnd=wproc.orgensemblesat('list.txt',11)
      you can see the information available with fdatacontrol.keys()
      where fdatacontrol['ohs'].values is the array with satellite observations of Hs
      and fdatacontrol['mhs'].values is the array with model results of Hs 
      for the control member.
     You can directly compose the whole ensemble (including control and
      perturbed members) with: np.c_[fdatacontrol['mhs'].values,fhs] where
      the shape[1] dimension contains the ensemble member.
      The ensemble mean will be np.nanmean(np.c_[fdatacontrol['mhs'].values,fhs],axis=1)
     If you want to save fdatacontrol,fhs,fwnd for latter use, you can do using pickle format:
     import pickle
     with open("WW3Altimeter_Data.pickle", 'wb') as file:
      pickle.dump(fdata, file)
      pickle.dump(fmhs, file)
      pickle.dump(fmwnd, file)
     and read later on using:
     with open('WW3Altimeter_Data.pickle', 'rb') as file:
      fdata = pickle.load(file)
      fmhs = pickle.load(file)
      fmwnd = pickle.load(file)
    '''

    if len(args) < 2:
        sys.exit(' Need two inputs with list name and number of members.')
    elif len(args) == 2:
        flist=args[0]; nmb=int(args[1])
    elif len(args) > 2:
        sys.exit(' Too many inputs')

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
                    'lat': np.array(f.variables['latitude'][:]).astype('double'),
                    'lon': np.array(f.variables['longitude'][:]).astype('double'),
                    'sid': np.array(f.variables['satelliteID'][:]).astype('double'),                    
                    'ohs': np.array(f.variables['obs_hs'][:]).astype('double')},
                    columns=['time','cycle','lat','lon','sid','ohs'])

                dfa = dfa[~dfa.duplicated(subset=['time','cycle','lat','lon','sid','ohs'])]
                dfa = dfa.reset_index(drop=True)

            dfb = pd.DataFrame({'time': np.array(f.variables['time'][:]).astype('double'),
                'cycle': np.array(f.variables['cycle'][:]).astype('double'),
                'lat': np.array(f.variables['latitude'][:]).astype('double'),
                'lon': np.array(f.variables['longitude'][:]).astype('double'),
                'sid': np.array(f.variables['satelliteID'][:]).astype('double'),
                'ohs': np.array(f.variables['obs_hs'][:]).astype('double')},
                columns=['time','cycle','lat','lon','sid','ohs'])

            dfb = dfb[~dfb.duplicated(subset=['time','cycle','lat','lon','sid','ohs'])]
            dfb = dfb.reset_index(drop=True)

            dfa = pd.merge(dfa, dfb)
            dfa = dfa[~dfa.duplicated(subset=['time','cycle','lat','lon','sid','ohs'])]
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
                'lat': np.array(f.variables['latitude'][:]).astype('double'),
                'lon': np.array(f.variables['longitude'][:]).astype('double'),
                'sid': np.array(f.variables['satelliteID'][:]).astype('double'),
                'ohs': np.array(f.variables['obs_hs'][:]).astype('double')},
                columns=['time','cycle','lat','lon','sid','ohs'])

            dfb = dfb[~dfb.duplicated(subset=['time','cycle','lat','lon','sid','ohs'])]
            # fast way to find the elements of dfb contained in dfa
            dfb = pd.merge(dfb, dfa)
            ind = np.array(dfb.index).astype('int')
            del dfb

            if j==0:
                df = pd.DataFrame({'time': np.array(f.variables['time'][:]).astype('double'),
                    'cycle': np.array(f.variables['cycle'][:]).astype('double'),
                    'lat': np.array(f.variables['latitude'][:]),
                    'lon': np.array(f.variables['longitude'][:]),
                    'ohs': np.array(f.variables['obs_hs'][:]),
                    'ownd': np.array(f.variables['obs_wnd'][:]),
                    'mhs': np.array(f.variables['model_hs'][:]),
                    'mwnd': np.array(f.variables['model_wnd'][:]),
                    'month': np.array(f.variables['month'][:]),
                    'distcoast': np.array(f.variables['distcoast'][:]),
                    'depth': np.array(f.variables['depth'][:]),
                    'GlobalOceansSeas': np.array(f.variables['GlobalOceansSeas'][:]),
                    'HighSeasMarineZones': np.array(f.variables['HighSeasMarineZones'][:]),
                    'cyclone': np.array(f.variables['cyclone'][:])},
                    columns=['time','cycle','lat','lon','ohs','ownd','mhs','mwnd','month','distcoast','depth','GlobalOceansSeas','HighSeasMarineZones','cyclone'])

                df = df.loc[ind]
                if c==0:
                    ocnames=np.array(f.variables['names_GlobalOceansSeas'][:]).astype('str')
                    hsmznames=np.array(f.variables['names_HighSeasMarineZones'][:]).astype('str')
                    cinfo=np.array(f.variables['cycloneinfo'][:]).astype('str')

            elif j==1:
                amhs=np.atleast_2d(np.array(f.variables['model_hs'][:])[ind])
                amwnd=np.atleast_2d(np.array(f.variables['model_wnd'][:])[ind])
            else:
                amhs=np.append(amhs,np.atleast_2d(np.array(f.variables['model_hs'][:])[ind]),axis=0)
                amwnd=np.append(amwnd,np.atleast_2d(np.array(f.variables['model_wnd'][:])[ind]),axis=0)

            f.close(); del ind,f
            c=c+1

        del dfa
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
        print('  Ok  '+"_".join(flist[i].split('_')[0:3]).split('/')[-1]+'  for all ensemble members')

    return fdata,fmhs,fmwnd
    print(" Done. wproc.orgensemblesat concluded.")
    # ------------------------

# Under development: spectral interpolation function, quality control, wind speed conversion to 10m



