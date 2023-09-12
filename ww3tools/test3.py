#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import netCDF4 as nc
import datetime
import time
from calendar import timegm
import sys
import os
import warnings
warnings.filterwarnings("ignore")

# Function to process data for a specific satellite
def process_satellite_data(dirs, sdname, sname, start_date, end_date):
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
                # Corrected file path using the updated satellite name and directory
                file_path = f"{dirs}/{sdname}/IMOS_SRS-Surface-Waves_MW_{sname}_FV02_{str(np.abs(j)).zfill(3)}{hem}-{str(k).zfill(3)}E-DM00.nc"
                fu = nc.Dataset(file_path)
            except FileNotFoundError:
                print(f"File {file_path} does not exist")
                continue

            st_strings = fu.variables['TIME'][:]
            st_strings = st_strings * 24. * 3600. + float(timegm(time.strptime('1985010100', '%Y%m%d%H')))
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
                    ast[ii:ii + len(valid_indices)] = np.array(st[valid_indices]).astype('float')
                    aslat[ii:ii + len(valid_indices)] = np.array(slat[valid_indices]).astype('float')
                    aslon[ii:ii + len(valid_indices)] = np.array(slon[valid_indices]).astype('float')
                    ahsk[ii:ii + len(valid_indices)] = np.array(hsk[valid_indices]).astype('float')
                    ahskcal[ii:ii + len(valid_indices)] = np.array(hskcal[valid_indices]).astype('float')
                    awnd[ii:ii + len(valid_indices)] = np.array(wnd[valid_indices]).astype('float')
                    awndcal[ii:ii + len(valid_indices)] = np.array(wndcal[valid_indices]).astype('float')
                    asig0knstd[ii:ii + len(valid_indices)] = np.array(sig0knstd[valid_indices]).astype('float')
                    aswhknobs[ii:ii + len(valid_indices)] = np.array(swhknobs[valid_indices]).astype('float')
                    aswhknstd[ii:ii + len(valid_indices)] = np.array(swhknstd[valid_indices]).astype('float')
                    aswhkqc[ii:ii + len(valid_indices)] = np.array(swhkqc[valid_indices]).astype('float')
                    ii = ii + len(valid_indices)
                else:
                    sys.exit('Small array to allocate the satellite data! Increase the power of the initial array (pia)')

                del st, slat, slon, hsk, hskcal, wnd, wndcal, sig0knstd, swhknobs, swhknstd, swhkqc
                fu.close()
                del fu

    print(f'Done reading and allocating satellite data for {sname}')

    # Check if any data was read and allocated
    if ii > 0:
        # Perform time averaging for every 10 minutes
        time_intervals = np.arange(start_date.timestamp(), end_date.timestamp(), 600)  # 10 minutes in seconds
        time_averaged_data = {
            'TIME': [],
            'LATITUDE': [],
            'LONGITUDE': [],
            'SWH_KU': [],
            'SWH_KU_CAL': [],
            'WSPD': [],
            'WSPD_CAL': [],
            'SIG0_KU_std_dev': [],
            'SWH_KU_num_obs': [],
            'SWH_KU_std_dev': [],
            'SWH_KU_quality_control': [],
        }

        for interval_start in time_intervals:
            interval_end = interval_start + 600  # 10 minutes in seconds
            interval_indices = np.where(np.logical_and(ast[:ii] >= interval_start, ast[:ii] < interval_end))[0]

            if len(interval_indices) > 0:
                # Calculate time-averaged values for each variable within the time interval
                time_averaged_data['TIME'].append(np.mean(ast[interval_indices]))
                time_averaged_data['LATITUDE'].append(np.mean(aslat[interval_indices]))
                time_averaged_data['LONGITUDE'].append(np.mean(aslon[interval_indices]))
                time_averaged_data['SWH_KU'].append(np.mean(ahsk[interval_indices]))
                time_averaged_data['SWH_KU_CAL'].append(np.mean(ahskcal[interval_indices]))
                time_averaged_data['WSPD'].append(np.mean(awnd[interval_indices]))
                time_averaged_data['WSPD_CAL'].append(np.mean(awndcal[interval_indices]))
                time_averaged_data['SIG0_KU_std_dev'].append(np.mean(asig0knstd[interval_indices]))
                time_averaged_data['SWH_KU_num_obs'].append(np.mean(aswhknobs[interval_indices]))
                time_averaged_data['SWH_KU_std_dev'].append(np.mean(aswhknstd[interval_indices]))
                time_averaged_data['SWH_KU_quality_control'].append(np.mean(aswhkqc[interval_indices]))

        return time_averaged_data

# Define the main program
if __name__ == '__main__':
    # Specify the data directory and start/end date
    dirs = "/scratch2/NCEPDEV/marine/Matthew.Masarik/dat/sat/AODN/altimeter"
    start_date = datetime.datetime(2012, 1, 1, 0, 0, 0)
    end_date = datetime.datetime(2013, 1, 2, 0, 0, 0)

    # Specify the satellite name and directory
    sdname = 'HY2'
    sname = 'HY-2'

    # Specify the power of the initial array
    pia = 10

    # Process the data for the specified satellite
    time_averaged_data = process_satellite_data(dirs, sdname, sname, start_date, end_date)

    if time_averaged_data is not None:
        # Save the time-averaged data in a new NetCDF file
        output_file = f"{sdname}_time_averaged_data.nc"
        with nc.Dataset(output_file, 'w', format='NETCDF4') as ncfile:
            # Define dimensions
            ncfile.createDimension('TIME', None)

            # Create variables
            time_var = ncfile.createVariable('TIME', 'f8', ('TIME',))
            lat_var = ncfile.createVariable('LATITUDE', 'f4', ('TIME',))
            lon_var = ncfile.createVariable('LONGITUDE', 'f4', ('TIME',))
            swh_ku_var = ncfile.createVariable('SWH_KU', 'f4', ('TIME',))
            swh_ku_cal_var = ncfile.createVariable('SWH_KU_CAL', 'f4', ('TIME',))
            wspd_var = ncfile.createVariable('WSPD', 'f4', ('TIME',))
            wspd_cal_var = ncfile.createVariable('WSPD_CAL', 'f4', ('TIME',))
            sig0_ku_std_dev_var = ncfile.createVariable('SIG0_KU_std_dev', 'f4', ('TIME',))
            swh_ku_num_obs_var = ncfile.createVariable('SWH_KU_num_obs', 'f4', ('TIME',))
            swh_ku_std_dev_var = ncfile.createVariable('SWH_KU_std_dev', 'f4', ('TIME',))
            swh_ku_qc_var = ncfile.createVariable('SWH_KU_quality_control', 'f4', ('TIME',))

            # Fill variables with time-averaged data
            time_var[:] = time_averaged_data['TIME']
            lat_var[:] = time_averaged_data['LATITUDE']
            lon_var[:] = time_averaged_data['LONGITUDE']
            swh_ku_var[:] = time_averaged_data['SWH_KU']
            swh_ku_cal_var[:] = time_averaged_data['SWH_KU_CAL']
            wspd_var[:] = time_averaged_data['WSPD']
            wspd_cal_var[:] = time_averaged_data['WSPD_CAL']
            sig0_ku_std_dev_var[:] = time_averaged_data['SIG0_KU_std_dev']
            swh_ku_num_obs_var[:] = time_averaged_data['SWH_KU_num_obs']
            swh_ku_std_dev_var[:] = time_averaged_data['SWH_KU_std_dev']
            swh_ku_qc_var[:] = time_averaged_data['SWH_KU_quality_control']

        print(f"Time-averaged data for {sname} saved in {output_file}")

