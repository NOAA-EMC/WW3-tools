#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import netCDF4 as nc
import datetime
import time
from calendar import timegm
import sys
import os
import warnings; warnings.filterwarnings("ignore")

# Function to process data for a specific satellite
def process_satellite_data(sdname, sname):
    # Define a range of latitudes and longitudes
    auxlat = np.arange(-90, 91, 1)
    auxlon = np.arange(0, 361, 1)  # Adjust the longitude range as needed

    # Initialize variables
    ast = np.zeros((10**pia), 'f')
    aslat = np.zeros((10**pia), 'f')
    aslon = np.zeros((10**pia), 'f')
    ahsk = np.zeros((10**pia), 'f')
    ahskcal = np.zeros((10**pia), 'f')
    awnd = np.zeros((10**pia), 'f')
    awndcal = np.zeros((10**pia), 'f')
    asig0knstd = np.zeros((10**pia), 'f')
    aswhknobs = np.zeros((10**pia), 'f')
    aswhknstd = np.zeros((10**pia), 'f')
    aswhkqc = np.zeros((10**pia), 'f')
    ii = 0

    # Loop through data
    for j in auxlat:
        for k in auxlon:
            if j >= 0:
                hem = 'N'
            else:
                hem = 'S'
            try:
                file_path = f"{dirs}/{sdname}/IMOS_SRS-Surface-Waves_MW_{sname}_FV02_{str(np.abs(j)).zfill(3)}{hem}-{str(k).zfill(3)}E-DM00.nc"
                fu = nc.Dataset(file_path)
            except FileNotFoundError:
                print(f"File {file_path} does not exist")
                continue

            st_strings = fu.variables['TIME'][:]
            st_strings=st_strings*24.*3600.+float(timegm( time.strptime('1985010100', '%Y%m%d%H') ))
            st = np.array([float(ts) for ts in st_strings])
            time_within_range = np.logical_and(st >= start_date.timestamp(), st <= end_date.timestamp())

            if np.sum(time_within_range) > 10:
                slat = fu.variables['LATITUDE'][:]
                slon = fu.variables['LONGITUDE'][:]
                wnd = fu.variables['WSPD'][:]
                wndcal = fu.variables['WSPD_CAL'][:]
                try:
                    hsk = fu.variables['SWH_KU'][:]
                    hskcal = fu.variables['SWH_KU_CAL'][:]
                    sig0knstd = fu.variables['SIG0_KU_std_dev'][:]
                    swhknobs = fu.variables['SWH_KU_num_obs'][:]
                    swhknstd = fu.variables['SWH_KU_std_dev'][:]
                    swhkqc = fu.variables['SWH_KU_quality_control'][:]
                except KeyError:
                    print('Error reading KU, picking KA')
                    hsk = fu.variables['SWH_KA'][:]
                    hskcal = fu.variables['SWH_KA_CAL'][:]
                    sig0knstd = fu.variables['SIG0_KA_std_dev'][:]
                    swhknobs = fu.variables['SWH_KA_num_obs'][:]
                    swhknstd = fu.variables['SWH_KA_std_dev'][:]
                    swhkqc = fu.variables['SWH_KA_quality_control'][:]

                valid_indices = np.where(time_within_range)[0]

                if ii + len(valid_indices) <= ast.shape[0]:
                    ast[ii:ii+len(valid_indices)] = np.array(st[valid_indices]).astype('float')
                    aslat[ii:ii+len(valid_indices)] = np.array(slat[valid_indices]).astype('float')
                    aslon[ii:ii+len(valid_indices)] = np.array(slon[valid_indices]).astype('float')
                    ahsk[ii:ii+len(valid_indices)] = np.array(hsk[valid_indices]).astype('float')
                    ahskcal[ii:ii+len(valid_indices)] = np.array(hskcal[valid_indices]).astype('float')
                    awnd[ii:ii+len(valid_indices)] = np.array(wnd[valid_indices]).astype('float')
                    awndcal[ii:ii+len(valid_indices)] = np.array(wndcal[valid_indices]).astype('float')
                    asig0knstd[ii:ii+len(valid_indices)] = np.array(sig0knstd[valid_indices]).astype('float')
                    aswhknobs[ii:ii+len(valid_indices)] = np.array(swhknobs[valid_indices]).astype('float')
                    aswhknstd[ii:ii+len(valid_indices)] = np.array(swhknstd[valid_indices]).astype('float')
                    aswhkqc[ii:ii+len(valid_indices)] = np.array(swhkqc[valid_indices]).astype('float')
                    ii = ii + len(valid_indices)
                else:
                    sys.exit('Small array to allocate the satellite data! Increase the power of the initial array (pia)')

                del st, slat, slon, hsk, hskcal, wnd, wndcal, sig0knstd, swhknobs, swhknstd, swhkqc
                fu.close()
                del fu

    print(f'Done reading and allocating satellite data for {sname}')
    
    # Check if any data was read and allocate
    if ii > 0:
        # Save netcdf
        ncfile = nc.Dataset('AltimeterData_' + sdname + '.nc', "w", format=fnetcdf)
        ncfile.createDimension('TIME', None)
        ncfile.createVariable('TIME', 'f4', ('TIME',))
        ncfile.variables['TIME'].units = 'days since 1985-01-01 00:00:00 UTC'
        ncfile.createVariable('LATITUDE', 'f4', ('TIME',))
        ncfile.createVariable('LONGITUDE', 'f4', ('TIME',))
        ncfile.createVariable('SWH_KU', 'f4', ('TIME',))
        ncfile.createVariable('SWH_KU_CAL', 'f4', ('TIME',))
        ncfile.createVariable('WSPD', 'f4', ('TIME',))
        ncfile.createVariable('WSPD_CAL', 'f4', ('TIME',))
        ncfile.createVariable('SIG0_KU_std_dev', 'f4', ('TIME',))
        ncfile.createVariable('SWH_KU_num_obs', 'f4', ('TIME',))
        ncfile.createVariable('SWH_KU_std_dev', 'f4', ('TIME',))
        ncfile.createVariable('SWH_KU_quality_control', 'f4', ('TIME',))

        ncfile.variables['TIME'][:] = ast[:ii]
        ncfile.variables['LATITUDE'][:] = aslat[:ii]
        ncfile.variables['LONGITUDE'][:] = aslon[:ii]
        ncfile.variables['SWH_KU'][:] = ahsk[:ii]
        ncfile.variables['SWH_KU_CAL'][:] = ahskcal[:ii]
        ncfile.variables['WSPD'][:] = awnd[:ii]
        ncfile.variables['WSPD_CAL'][:] = awndcal[:ii]
        ncfile.variables['SIG0_KU_std_dev'][:] = asig0knstd[:ii]
        ncfile.variables['SWH_KU_num_obs'][:] = aswhknobs[:ii]
        ncfile.variables['SWH_KU_std_dev'][:] = aswhknstd[:ii]
        ncfile.variables['SWH_KU_quality_control'][:] = aswhkqc[:ii]
        ncfile.close()
        print('Data saved to AltimeterData_' + sdname + '.nc')
    else:
        print('No valid data found for ' + sdname)


if __name__ == "__main__":

    # Define constants and parameters
    fnetcdf = "NETCDF4"
    pia = 10
    dirs = '/scratch2/NCEPDEV/marine/Matthew.Masarik/dat/sat/AODN/altimeter'

    # Satellite missions available at AODN dataset
    sdnames = np.array(['SENTINEL3A'])
    snames = np.array(['SENTINEL-3A'])

    # Define the time range
    start_date = datetime.datetime(2017, 12, 1, 0, 0)  # December 1, 2015, at midnight
    end_date = datetime.datetime(2018, 1, 31, 23, 59)  # January 31, 2016, just before midnight

    # Process data for each satellite
    for s in range(len(sdnames)):
        print(f"Processing data for satellite {snames[s]} [{s + 1}/{len(sdnames)}]")
        process_satellite_data(sdnames[s], snames[s])


