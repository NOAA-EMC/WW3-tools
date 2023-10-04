#Description:

#The provided Python script (SatGlobal_Altimeter_new.py) is designed to process satellite data based on the input parameters you provide.
#The batch script sets these parameters dynamically and submits the job for processing. Make sure to modify the batch script
#to match your specific requirements, such as the satellite name, time-averaging settings, and paths to input and output directories.

#In the batch script one should define:
#Set your desired start and end dates using the start_date and end_date variables.
#Specify the input directory (input_directory) where your satellite data files are located.
#Specify the output directory (output_directory) where you want to save the processed data.
#Set the satellite name (satellite_name) and other time-averaging parameters (time_averaging, averaging_time) as needed.


#DEPENDENCIES:
# Before using this script, make sure you have the necessary Python libraries installed:
# - numpy
# - netCDF4
# - datetime
# - time
# - calendar
# - sys
# - os
#AODN altimeter data previously downloaded (see wfetchsatellite_AODN_Altimeter.sh)

#Author:
#Maryam Mohammadpour, Ricardo Campos

#Contact:
#Ghazal.Mohammadpour@noaa.gov

#---------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------

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
import argparse

warnings.filterwarnings("ignore")

# Function to process data for a specific satellite
def process_satellite_data(dirs, sdname, sname, start_date, end_date, enable_time_averaging, time_averaging_duration):
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
        if enable_time_averaging:
            # Perform time averaging for every specified duration
            time_intervals = np.arange(start_date.timestamp(), end_date.timestamp(), time_averaging_duration)
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
                interval_end = interval_start + time_averaging_duration
                interval_indices = np.where(np.logical_and(ast[:ii] >= interval_start, ast[:ii] < interval_end))[0]

                if len(interval_indices) > 0:
                    # Calculate time-averaged values for each variable within the time interval
                    time_averaged_data['TIME'].append(np.nanmean(ast[interval_indices]))
                    time_averaged_data['LATITUDE'].append(np.nanmean(aslat[interval_indices]))
                    time_averaged_data['LONGITUDE'].append(np.nanmean(aslon[interval_indices]))
                    time_averaged_data['SWH_KU'].append(np.nanmean(ahsk[interval_indices]))
                    time_averaged_data['SWH_KU_CAL'].append(np.nanmean(ahskcal[interval_indices]))
                    time_averaged_data['WSPD'].append(np.nanmean(awnd[interval_indices]))
                    time_averaged_data['WSPD_CAL'].append(np.nanmean(awndcal[interval_indices]))
                    time_averaged_data['SIG0_KU_std_dev'].append(np.nanmean(asig0knstd[interval_indices]))
                    time_averaged_data['SWH_KU_num_obs'].append(np.nanmean(aswhknobs[interval_indices]))
                    time_averaged_data['SWH_KU_std_dev'].append(np.nanmean(aswhknstd[interval_indices]))
                    time_averaged_data['SWH_KU_quality_control'].append(np.nanmean(aswhkqc[interval_indices]))

            # Convert lists to numpy arrays
            for key, value in time_averaged_data.items():
                time_averaged_data[key] = np.array(value)

            return time_averaged_data

        else:
            # If not performing time averaging, return the raw data
            raw_data = {
                'TIME': ast[:ii],
                'LATITUDE': aslat[:ii],
                'LONGITUDE': aslon[:ii],
                'SWH_KU': ahsk[:ii],
                'SWH_KU_CAL': ahskcal[:ii],
                'WSPD': awnd[:ii],
                'WSPD_CAL': awndcal[:ii],
                'SIG0_KU_std_dev': asig0knstd[:ii],
                'SWH_KU_num_obs': aswhknobs[:ii],
                'SWH_KU_std_dev': aswhknstd[:ii],
                'SWH_KU_quality_control': aswhkqc[:ii],
            }
            return raw_data
    else:
        print(f'No valid data found for {sname}')
        return None

# Function to parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser(description='Process satellite data')
    parser.add_argument('--time-averaging', action='store_true', help='Enable time averaging')
    parser.add_argument('--duration', type=int, default=600, help='Time averaging duration in seconds')
    parser.add_argument('--start-date', required=True, help='Start date in the format "YYYY-MM-DD HH:MM:SS"')
    parser.add_argument('--end-date', required=True, help='End date in the format "YYYY-MM-DD HH:MM:SS"')
    parser.add_argument('-i', '--input', required=True, help='Input directory path')
    parser.add_argument('-o', '--output', default='output.nc', help='Output file path')
    parser.add_argument('--sname', required=True, help='Satellite name')
    parser.add_argument('--sdname', required=True, help='Satellite directory name')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()

    # Parse start and end dates
    start_date = datetime.datetime.strptime(args.start_date, '%Y-%m-%d %H:%M:%S')
    end_date = datetime.datetime.strptime(args.end_date, '%Y-%m-%d %H:%M:%S')

    # Specify the data directory
    dirs = args.input  # Use the input path provided as an argument

    # Specify the satellite name and directory
    sdname = args.sdname  # Use the satellite directory provided as an argument
    sname = args.sname  # Use the satellite name provided as an argument

    # Specify the power of the initial array
    pia = 10

    # Process the data for the specified satellite with time averaging as specified by the user
    data = process_satellite_data(dirs, sdname, sname, start_date, end_date, args.time_averaging, args.duration)

    if data is not None:
        if args.time_averaging:
            output_file = f"{sdname}_time_averaged_data_{start_date.strftime('%Y%m%d%H%M%S')}_{end_date.strftime('%Y%m%d%H%M%S')}_{args.duration}s.nc"
        else:
            output_file = f"{sdname}_raw_data.nc"

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

            # Fill variables with data as before
            time_var[:] = data['TIME']
            lat_var[:] = data['LATITUDE']
            lon_var[:] = data['LONGITUDE']
            swh_ku_var[:] = data['SWH_KU']
            swh_ku_cal_var[:] = data['SWH_KU_CAL']
            wspd_var[:] = data['WSPD']
            wspd_cal_var[:] = data['WSPD_CAL']
            sig0_ku_std_dev_var[:] = data['SIG0_KU_std_dev']
            swh_ku_num_obs_var[:] = data['SWH_KU_num_obs']
            swh_ku_std_dev_var[:] = data['SWH_KU_std_dev']
            swh_ku_qc_var[:] = data['SWH_KU_quality_control']

            # Set attributes for time average start time, end time, and duration
            ncfile.setncattr('TimeAverageStartTime', str(start_date))
            ncfile.setncattr('TimeAverageEndTime', str(end_date))
            ncfile.setncattr('TimeAverageDuration', str(args.duration))

        print(f"Data for {sname} saved in {output_file}")


