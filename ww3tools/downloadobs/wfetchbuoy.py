#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
wfetchbuoy.py

VERSION AND LAST UPDATE:
 v1.1  01/08/2021
 v1.2  04/04/2022
 v1.3  01/07/2022
 v1.4  01/17/2023
 v1.5  01/18/2023
 v1.6  01/25/2023

PURPOSE:
 This module downloads NDBC metocean data. One file per buoy.

 NDBC data can be downloaded for a user-provided time period (files are divided by year) 
 Information about NDBC formats, netcdf and stdmet, can be found at:
  https://dods.ndbc.noaa.gov/thredds/fileServer/data/
  https://www.ndbc.noaa.gov/rsa.shtml
  https://www.ndbc.noaa.gov/measdes.shtml
  https://www.ndbc.noaa.gov/stndesc.shtml
 The NDBC stdmet function (ndbc_stdmet) retrieves only integrated wave parameters,
  while the function ndbc_nc downloads both wave-parameters and spectra.
 NOAA National Data Buoy Center:
  https://www.ndbc.noaa.gov/

 See also https://pypi.org/project/ndbc-api/

 For Copernicus buoy data:
  https://help.marine.copernicus.eu/en/articles/7949409-copernicus-marine-toolbox-introduction
  https://data.marine.copernicus.eu/product/INSITU_GLO_WAV_DISCRETE_MY_013_045/files?subdataset=cmems_obs-ins_glo_wav_my_na_PT1H_202311--ext--history

USAGE:
 functions
   ndbc_nc : downloads NDBC buoy dataset integrated parameters and spectrum, in netcdf format
   ndbc_stdmet : downloads time-series of NDBC buoy dataset of integrated parameters in text stdmet format

 A list of buoy IDs is required and final output paths can be given.
 For ndbc_nc and ndbc_stdmet, two inputs of initial and final years must be entered.
 Explanation for each function is contained in the headers.

OUTPUT:
 Formatted and saved files for all available buoy stations given in a list.
 One file per buoy.

DEPENDENCIES:
 See setup.py and the imports below.

AUTHOR and DATE:
 01/07/2021: Ricardo M. Campos, first version considering NDBC stdmet format only.
 01/08/2021: NDBC netcdf format, Jake Campbell and Ali Abdolali, first version (originally retrieve_ndbc_nc.py).
 04/04/2022: stdmet format function updates, Ricardo M. Campos.
 01/07/2022: Ricardo Campos, included wget option in ndbc_nc and it reads a list of buoys.
 01/17/2023: Ricardo Campos, code updated to fit into ww3_tools.
 01/18/2023: Ricardo Campos, NDBC nc and stdmet formats adapted to functions and unified in the same module.
 01/25/2023: Ricardo Campos, NDBC and Copernicus funtions adapted and included in the same module, 
  i.e., get_buoydata_copernicus.py and get_ndbc.py became wfetchbuoy.py.

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import os
import sys
import numpy as np
import netCDF4
import pandas as pd
from pylab import *
import io
import urllib
from urllib import request
from pathlib2 import Path
import xarray as xr
import time
from tqdm import *
import warnings; warnings.filterwarnings("ignore")

#  ----- NDBC BUOY DATA -----

def ndbc_nc(*args):
    """
    wfetchbuoy.ndbc_nc , adapted from retrieve_ndbc_nc.py written by Jake Campbell and Ali Abdolali.

    PURPOSE:
     This function downloads the historical NDBC data (both wave parameters and spectra)
      for a user provided time period (initial and final year)
     NDBC buoy data in netcdf format is obtained:
      https://dods.ndbc.noaa.gov/thredds/fileServer/data/
     for both wave-parameters and spectra:
     https://dods.ndbc.noaa.gov/thredds/fileServer/data/stdmet/
     https://dods.ndbc.noaa.gov/thredds/fileServer/data/swden/
     The download is primarily done with urllib and xarray but, when it fails,
      wget is then used (linux system dependent).
     The url used for download, when entering a past time interval, refers to 
     quality-controlled observations of historical data (archive). 
     Very recent data (near-real time) can be obtained by entering 9999 for the 
      initial and final years. In this case, the last 10 years will be downloaded, up to 
      the last hour (if available).
     NOAA National Data Buoy Center:
      https://www.ndbc.noaa.gov/

    USAGE:
     Five arguments required (first three mandatory):
      first year (or 9999 for the last 10 years covering up to near-real-time data)
      last year (or 9999 for the last 10 years covering up to near-real-time data)
      station list can be an array/list or a text file (full path) with one column (.txt) or a ww3 point input format file (.dat)
      output dir wave parameters
      output dir wave spectra
     If output dirs are not provided, data will be saved in the local working directory.
     Examples:
      from ww3tools.downloadobs import wfetchbuoy
      wfetchbuoy.ndbc_nc(2012,2022,'/home/name/reference_data/buoys/allbstations.txt','/home/name/reference_data/buoys/NDBC/wparam','/home/name/reference_data/buoys/NDBC/spec')
      wfetchbuoy.ndbc_nc(9999,9999,'atlanticbuoys.txt')
      wfetchbuoy.ndbc_nc(1990,2020,[41004,41047])

    OUTPUT:
     Formatted and saved netcdf files for all available NDBC buoy stations
      for both wave parameters and spectra.
     One file per buoy.

    """

    saveSpectra = os.getcwd()+'/' ; saveWave = os.getcwd()+'/'
    if len(args) < 3 :
        sys.exit(' At least three arguments (first year, last year, and buoy_list) must be informed.')
    if len(args) >= 3 :
        iyear = int(args[0])
        fyear = int(args[1])
        lstations=args[2]
    if len(args) >= 4:
        saveWave = str(args[3])
        if saveWave[-1]!='/':
            saveWave=saveWave+'/'

        saveSpectra=np.copy(saveWave)

    if len(args) >= 5:
        saveSpectra = str(args[4])
        if saveSpectra[-1]!='/':
            saveSpectra=saveSpectra+'/'

    if len(args) > 5:
        sys.exit(' Too many inputs')

    # Basing station names off current list of NDBC bouys, reading them from text file
    if type(lstations)==str:
        if lstations.split('.')[-1] == 'dat':
            # WAVEWATCHIII point input format
            dfabs = pd.read_csv(lstations, comment='$'); stationList=[]
            for i in range(0,dfabs.values.shape[0]-1):
                aux=str(dfabs.values[i]).split()
                stationList=np.append(stationList,str(aux[3][1::]).split("'")[0])
        else:
            # List (single column) of buoy IDs
            f = open(lstations, "r")
            stationList = f.read().splitlines()
            # List of lon lat buoyID
            # for i in range(0,len(stationList)):
            #    stationList[i]=stationList[i].split()[-1]

    else:
        stationList=np.array(lstations).astype('str')
    ###############################################################################

    # While loop to loop through Station List for both wave and spectra NetCDF 
    # files and if they exist for given year (tqdm is simply a progress bar)
    with tqdm(total=int(len(stationList)*((fyear-iyear)+1)), desc=" -Downloading NDBC Files (.nc)", position=0, leave=True) as pbar:

        for i in range(0,len(stationList)):

            for yr in range(iyear,fyear+1):
                start = str(yr)+'-01-01'
                end = str(yr+1)+'-01-01'

                time.sleep(.01)
                stat = stationList[i]
                webWave = 'https://dods.ndbc.noaa.gov/thredds/fileServer/data/stdmet/'+str(stat)+'/'+str(stat)+'h'+str(yr)+'.nc'
                webSpectra = 'https://dods.ndbc.noaa.gov/thredds/fileServer/data/swden/'+str(stat)+'/'+str(stat)+'w'+str(yr)+'.nc'
                print(str(stat)+' '+str(yr))
                # Set save location for each
                saveloc = saveSpectra + str(stat) +'w'+ str(yr) +'.nc'
                saveloc2 = saveWave + str(stat) +'h'+ str(yr) +'.nc'

                # Perform pull with check to see if spectra file exists. If it does, skip and look for next station.
                if os.path.isfile(saveloc):
                    pass         
                else:
                    try:
                        reqSpectra = request.Request(webSpectra)
                        respS = request.urlopen(reqSpectra)
                    except:
                        pass
                    else:            
                        try:
                            ds_s = xr.open_dataset(io.BytesIO(respS.read()))
                        except :
                            # Linux system only
                            os.system('wget --no-check-certificate --no-proxy -l1 -H -t1 -nd -N -np -erobots=off --tries=3 '+webSpectra)
                            os.system('mv -f '+str(stat)+'w'+str(yr)+'.nc '+saveSpectra)
                            print('Warining with '+webSpectra+'  downloaded with wget ...')
                        else:
                            # Filter for specific date range before writing to NetCDF to save disk space.
                            if yr != 9999:
                                ds_s = ds_s.sel(time=slice(start, end)) 

                            # Check to see if dataset is empty, if so, do not save it.
                            ds_s_np = np.array(ds_s['time'])
                            if ds_s_np.shape[0] == 0:
                                pass
                            else:
                                ds_s.to_netcdf(saveloc)
                
                # Perform pull with check to see if wave file exists. If it does, skip and look for next station.    
                if os.path.isfile(saveloc2):
                    pass
                else:
                    try:
                        reqWave = request.Request(webWave)
                        respW = request.urlopen(reqWave)
                    except:
                        pass
                    else:
                        try:
                            ds_w = xr.open_dataset(io.BytesIO(respW.read()))
                        except :
                            os.system('wget --no-check-certificate --no-proxy -l1 -H -t1 -nd -N -np -erobots=off --tries=3 '+webWave)
                            os.system('mv -f '+str(stat)+'h'+str(yr)+'.nc '+saveWave)
                            print('Warining with '+webWave+'  downloaded with wget')
                        else:
                            # Can filter for specific days before writing to NetCDF.
                            if yr != 9999:
                                ds_w = ds_w.sel(time=slice(start, end))

                            # Check to see if dataset is empty, if so, do not save it.
                            ds_w_np = np.array(ds_w['time'])
                            if ds_w_np.shape[0] == 0:
                                pass
                            else:
                                ds_w.to_netcdf(saveloc2)

                pbar.update(1)

        pbar.close()


def ndbc_stdmet(*args):
    """
    wfetchbuoy.ndbc_stdmet

    PURPOSE:
     Download the historical NDBC buoy data in text stdmet format (time-series of wave parameters only)
      https://www.ndbc.noaa.gov/rsa.shtml
      https://www.ndbc.noaa.gov/measdes.shtml
      https://www.ndbc.noaa.gov/stndesc.shtml
     The url used for download refers to quality-controlled observations,
      and the code is meant for historical data (archive). Recent data (last
      week or month) is not included.
     For near-real time data, a different url source and script must be used.
     NOAA National Data Buoy Center:
      https://www.ndbc.noaa.gov/

    USAGE:
     Four arguments required (first three mandatory): firstYear, lastYear, station list, output dir
     station list can be an array/list or a text file (full path) with one column (.txt) or a ww3 point input format file (.dat)
     Examples:
      from ww3tools.downloadobs import wfetchbuoy
      wfetchbuoy.ndbc_stdmet(2010,2020,'allbstations.txt')
      wfetchbuoy.ndbc_stdmet(2012,2022,'/home/name/reference_data/buoys/allbstations.txt','/home/name/reference_data/buoys/NDBC')
      wfetchbuoy.ndbc_stdmet(1990,2020,[41004,41047])

    OUTPUT:
     Text files with NDBC metocean data for the buoys listed in 
      the station list (third argument).
     One file per buoy.

    """

    odir = os.getcwd()+'/'
    if len(args) < 3 :
        sys.exit(' At least three arguments (first year, last year, and buoy_list) must be informed.')
    if len(args) >= 3 :
        iyear = str(args[0])
        fyear = str(args[1])
        lstations=args[2]
    if len(args) >= 4:
        odir = str(args[3])
        if odir[-1]!='/':
            odir=odir+'/'

    if len(args) > 4:
        sys.exit(' Too many inputs')

    # Read list of buoys
    #lndbc = open('allbstations.txt'); content = lndbc.readlines()
    if type(lstations)==str:
        if lstations.split('.')[-1] == 'dat':
            # WAVEWATCHIII point input format
            dfabs = pd.read_csv(lstations, comment='$'); stationList=[]
            for i in range(0,dfabs.values.shape[0]-1):
                aux=str(dfabs.values[i]).split()
                stationList=np.append(stationList,str(aux[3][1::]).split("'")[0])
        else:
            # List (single column) of buoy IDs
            f = open(lstations, "r")
            stationList = f.read().splitlines()
            # List of lon lat buoyID
            # for i in range(0,len(stationList)):
            #    stationList[i]=stationList[i].split()[-1]
    else:
        stationList=np.array(lstations).astype('str')
    ###############################################################################

    with tqdm(total=int(len(stationList)*((int(fyear)-int(iyear))+1)), desc=" -Downloading NDBC Files (stdmet)", position=0, leave=True) as pbar:

        # Dowload each buoy/station for each year
        for st in range(0,np.size(stationList)):

            namest = str(stationList[st]).split()[0].lower()
            # tb=0
            for year in range(int(iyear),int(fyear)+1):

                url = 'http://www.ndbc.noaa.gov/view_text_file.php?filename='+namest+'h'+repr(year)+'.txt.gz&dir=data/historical/stdmet/'

                try:
                    response = request.urlopen(url) 
                except :
                    print(url+"   does not exist"); print(" ")
                else:

                    data = response.read() 
                    text = data.decode('utf-8')
                    del data, response, url

                    if len(text.splitlines()) >10:
                        text_file = open(odir+"/NDBC_stdmet_"+namest.upper()+"h"+repr(year)+".txt", "a")
                        text = text.replace("#YY","YY")
                        text_file.write(text)
                        text_file.close()
                        # tb=tb+1

                    del text
                    print(" NDBC buoy "+namest.upper()+" "+repr(year)+" (stdmet format) successfully downloaded."); print(" ")
                    
                pbar.update(1)

        del namest
        pbar.close()

