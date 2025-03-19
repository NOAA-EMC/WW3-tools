"""

# This Python script is designed for processing and interpolating oceanographic model data and satellite
# time-averaged observations, focusing on significant wave height (HTSGW, swh) and wind speed (WIND, ws) parameters.
# It supports handling data in both GRIB2 and NetCDF formats, aligning model outputs with satellite observations
# to facilitate comparison and analysis for oceanographic research.

# The script performs several key operations:
# 1. Data Loading: It loads model data from a specified directory based on a naming pattern and satellite data from a given file.
# 2. Time Alignment: Converts model and satellite data times to UNIX time for consistent time referencing and aligns them.
# 3. Data Preparation and Interpolation: Applies an ocean mask to model data and interpolates it to match satellite observation points.
# 4. Data Merging: Combines interpolated model data with satellite observations into a single xarray Dataset.
# 5. Output: Saves the combined dataset as a NetCDF file for further analysis or comparison studies.

# For the NetCDF data interpolation, the process involves:
# - Reading the model data from NetCDF files according to a specified pattern.
# - Aligning model data times with satellite observation times by converting all times to a common UNIX format.
# - Applying a spatial mask to filter out land data points and focus on oceanic areas for interpolation.
# - Using a RegularGridInterpolator for spatially interpolating model data (wind speed and significant wave height) to the satellite observation points.
# - Creating a comprehensive dataset that merges interpolated model data with the original satellite observations, ensuring that all data shares a common temporal and spatial frame of reference.
# - Saving this merged dataset as a new NetCDF file, ready for analysis, visualization, or further study.

# Requirements:
# - xarray: For handling data arrays and NetCDF files.
# - numpy: For numerical operations.
# - glob: For matching file path patterns.
# - re: For regular expression operations.
# - scipy.interpolate: For performing data interpolation.
# - pandas: For data manipulation.
# - cfgrib: For reading GRIB files.

# This script is intended for use by oceanographic data analysts and researchers who need to compare model output with satellite measurements.

import sys
import glob
import numpy as np
import xarray as xr
import re
from scipy.interpolate import RegularGridInterpolator
import cfgrib
import os
import pandas as pd

# Function to interpolate GRIB2 format model data
def interpolate_grib2(grib_data_directory, grib_data_pattern, satellite_file, output_file, satellite_name, ic_time):

    Interpolates GRIB2 format model data to satellite observation points and merges them into a NetCDF file.

    Parameters:
    - grib_data_directory (str): Directory containing GRIB2 files.
    - grib_data_pattern (str): Pattern to match GRIB2 files.
    - satellite_file (str): Satellite data file path.
    - output_file (str): Output NetCDF file path.
    - satellite_name (str): Name of the satellite for data labeling.


    # Implementation details here...
    pass

# Function to interpolate NetCDF format model data
def interpolate_netcdf(model_data_directory, model_data_pattern, satellite_file, output_file):

    Interpolates NetCDF format model data to satellite observation points and merges them into a NetCDF file.

    Parameters:
    - model_data_directory (str): Directory containing NetCDF files.
    - model_data_pattern (str): Pattern to match NetCDF files.
    - satellite_file (str): Satellite data file path.
    - output_file (str): Output NetCDF file path.

    This function performs the following steps:
    - Reads and concatenates model data from multiple NetCDF files based on the provided pattern.
    - Converts time data to UNIX format for both model and satellite datasets to ensure time alignment.
    - Applies a spatial mask to select oceanic areas, filtering out data points over land.
    - Interpolates model data to the geographic and temporal points of satellite observations using a RegularGridInterpolator.
    - Merges the interpolated model data with satellite observations, ensuring consistency in time and space.
    - Saves the merged dataset as a new NetCDF file, providing a comprehensive view of model predictions and actual satellite measurements.

    # Implementation details here...
    pass



usage: ProcSat_interpolation.py [-h] -t TYPEFILE -d DATADIR -p PATTERN -s SATFILE -o OUTDIR -f FILEOUT -m MODEL

options:
  -h, --help            show this help message and exit
  -t TYPEFILE, --typefile TYPEFILE
                        Type of Model File, 'grib2' or 'nc'
  -d DATADIR, --datadir DATADIR
                        Data Directory for Model Files
  -p PATTERN, --pattern PATTERN
                        Pattern of Model Files
  -s SATFILE, --satfile SATFILE
                        Satellite File
  -o OUTDIR, --outdir OUTDIR
                        Directory Path for Output
  -f FILEOUT, --fileout FILEOUT
                        Name of Output File
  -m MODEL, --model MODEL
                        String Identifier of Model

# Author: Ghazal Mohammadpour
# Email: ghazal.mohammadpour@noaa.gov
# This script modification ensures alignment with the provided requirements and specifications for processing oceanographic model and satellite data.


"""

import sys
import glob
import numpy as np
import xarray as xr
import re
from scipy.interpolate import RegularGridInterpolator
import cfgrib
import os
import pandas as pd
import argparse
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)


# Function to process_grib_files from a GRIB file
def process_grib_files(file):  
    ds = cfgrib.open_dataset(file, backend_kwargs={'indexpath': ''})
    if ds is None:
        raise ValueError(f"No dataset found in the file {file}")
    time = ds.coords['time'].values
    return time.item() if time.size == 1 else time

def extract_reference_times_grib(files):
    model_data_list = []
    valid_times_unix = []

    for file in files:
        ds = cfgrib.open_datasets(file, backend_kwargs={'indexpath': ''})
        ws = ds[1]['ws'] if 'ws' in ds[1] else xr.full_like(ds[1]['latitude'], np.nan, dtype=np.float32)
        swh = ds[1]['swh'] if 'swh' in ds[1] else xr.full_like(ds[1]['latitude'], np.nan, dtype=np.float32)
        valid_time_unix = (ds[1].valid_time.values - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')

        valid_times_unix.append(valid_time_unix.item())  # Append valid time

        data = xr.Dataset({
            'wind_speed': ws,
            'significant_wave_height': swh,
            'longitude': ds[1]['longitude'],
            'latitude': ds[1]['latitude'],
            'time': ds[1]['time'],
            'step': ds[1]['step']
        })

        model_data_list.append(data)

    return model_data_list, np.array(valid_times_unix)

def interpolate_grib2(data_directory, data_pattern, satellite_file, output_file, model_name):
    found_files = glob.glob(os.path.join(data_directory, data_pattern))
    if not found_files:
        print("No files found at:", found_files)
        raise ValueError("No GRIB files found in the specified directory matching the pattern.")

    found_files.sort()
    filtered_grib_files = found_files

    model_data_list, valid_times_unix = extract_reference_times_grib(filtered_grib_files)
    model_data_combined = xr.concat(model_data_list, dim='time')

    # Added longitude conversion
    model_longitudes = model_data_combined.longitude.values
    if np.any(model_longitudes < 0):  # Check if any longitudes are negative
        model_longitudes[model_longitudes < 0] += 360
        model_data_combined = model_data_combined.assign_coords(longitude=(model_longitudes))

    satellite_data = xr.open_dataset(satellite_file)
    satellite_times_unix = np.array((satellite_data['time'].values - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')).astype('double')

    if 'sat_name' in satellite_data.variables:
        satellite_name = satellite_data['sat_name'].values[0]
    elif 'satellite_name' in satellite_data.attrs:
        satellite_name = satellite_data.attrs['satellite_name']
    else:
        satellite_name = "unknown_satellite"

    min_model_time = min(valid_times_unix)
    max_model_time = max(valid_times_unix)

    time_mask = (satellite_times_unix >= min_model_time) & (satellite_times_unix <= max_model_time)
    filtered_satellite_times_unix = satellite_times_unix[time_mask]

    hs = satellite_data['hs'].values[time_mask]
    wsp_cal = satellite_data['wsp_cal'].values[time_mask]
    satellite_longitudes = satellite_data['longitude'].values[time_mask]
    satellite_longitudes[satellite_longitudes < 0] += 360   #added
    satellite_latitudes = satellite_data['latitude'].values[time_mask]
    hs_cal = satellite_data['hs_cal'].values[time_mask]
    wsp = satellite_data['wsp'].values[time_mask]

    time_dataarray = xr.DataArray(filtered_satellite_times_unix, dims=['time'], name='time', attrs={
        'standard_name': 'time',
        'units': 'seconds since 1970-01-01 00:00:00',
        'calendar': 'standard',
        'axis': 'T'
    })

    interpolated_ws_values, model_times_ws = prepare_and_interpolate(model_data_combined, 'wind_speed', satellite_latitudes, satellite_longitudes, filtered_satellite_times_unix, valid_times_unix)
    interpolated_swh_values, model_times_swh = prepare_and_interpolate(model_data_combined, 'significant_wave_height', satellite_latitudes, satellite_longitudes, filtered_satellite_times_unix, valid_times_unix)

    interpolated_dataset = xr.Dataset({
        'time': time_dataarray,
        'latitude': xr.DataArray(satellite_latitudes, coords={'time': filtered_satellite_times_unix}, dims=['time'], name='latitude').assign_attrs(units='degree_north'),
        'longitude': xr.DataArray(satellite_longitudes, coords={'time': filtered_satellite_times_unix}, dims=['time'], name='longitude').assign_attrs(units='degree_east'),
        'model_hs': xr.DataArray(interpolated_swh_values, coords={'time': filtered_satellite_times_unix}, dims=['time'], name='model_hs').assign_attrs(units='m'),
        'model_wnd': xr.DataArray(interpolated_ws_values, coords={'time': filtered_satellite_times_unix}, dims=['time'], name='model_wnd').assign_attrs(units='m'),
        'obs_hs': xr.DataArray(hs, coords={'time': filtered_satellite_times_unix}, dims=['time'], name='obs_hs').assign_attrs(units='m'),
        'obs_hs_cal': xr.DataArray(hs_cal, coords={'time': filtered_satellite_times_unix}, dims=['time'], name='obs_hs_cal').assign_attrs(units='m'),
        'obs_wnd': xr.DataArray(wsp, coords={'time': filtered_satellite_times_unix}, dims=['time'], name='obs_wnd').assign_attrs(units='m/s'),
        'obs_wnd_cal': xr.DataArray(wsp_cal, coords={'time': filtered_satellite_times_unix}, dims=['time'], name='obs_wnd_cal').assign_attrs(units='m/s'),
        'fcst_hr': xr.DataArray(model_times_swh, dims=['time'], name='fcst_hr').assign_attrs(units='seconds'),
    })

    interpolated_dataset.attrs['satellite_name'] = satellite_name
    interpolated_dataset.attrs['model_name'] = model_name
    initial_condition_time = process_grib_files(filtered_grib_files[0])
    initial_condition_time_unix = int(initial_condition_time) // 10**9
    fcst_hr2 = np.full_like(interpolated_ws_values, np.nan)  # Initialize with NaN
    valid_indices = ~np.isnan(interpolated_ws_values)
    fcst_hr2[valid_indices] = (filtered_satellite_times_unix[valid_indices] - initial_condition_time_unix) / 3600  # Calculate only for valid ws


    interpolated_dataset['fcst_hr'] = xr.DataArray(fcst_hr2, dims=['time'], name='fcst_hr2').assign_attrs(description="Forecast hour relative to initial condition time", units='hours')
    interpolated_dataset.attrs['initial_condition_time'] = float(initial_condition_time_unix)
    interpolated_dataset.attrs['initial_condition_time_units'] = "seconds since 1970-01-01 00:00:00"  # Specify units separately

    interpolated_dataset.to_netcdf(output_file, format='NETCDF4')

def spatial_mask(satellite_latitudes, satellite_longitudes,model_data):
    min_lat_val = model_data['latitude'].min().values.item()
    max_lat_val = model_data['latitude'].max().values.item()
    min_lon_val = model_data['longitude'].min().values.item()
    max_lon_val = model_data['longitude'].max().values.item()

    within_lat_bounds = (satellite_latitudes >= min_lat_val) & (satellite_latitudes <= max_lat_val)
    within_lon_bounds = (satellite_longitudes >= min_lon_val) & (satellite_longitudes <= max_lon_val)
    return within_lat_bounds & within_lon_bounds


# NETCDF interpolator
def interpolate_netcdf(model_data_directory, model_data_pattern, satellite_file, output_file, model_name):
    # Using glob to find all files in the directory that match the pattern
    model_data_files = glob.glob(os.path.join(model_data_directory, model_data_pattern))
    # Check if any files are found
    if not model_data_files:
        raise ValueError("No files found in the specified directory matching the pattern.")

    # Sort the files to ensure they are in the correct order
    model_data_files.sort()

    def extract_reference_times_netcdf(file_path):
        ds = xr.open_dataset(file_path)
        reference_time = ds.time.attrs['reference_time']
        reference_date = ds.time.attrs['reference_date']
        ds.close()  # Close the dataset after extraction
        return reference_time, reference_date

    model_reference_time_unix, reference_time_unix = extract_reference_times_netcdf(model_data_files[0])

    # Load and concatenate the model data files
    model_data_list = [xr.open_dataset(file) for file in model_data_files]
    model_data_combined = xr.concat(model_data_list, dim='time')
    model_data_combined.close()  # Close the files after loading

    # Convert longitudes from -180 to 180 to 0 to 360, if necessary
    model_longitudes = model_data_combined['longitude'].values
    if model_longitudes.min() < 0:
        model_longitudes[model_longitudes < 0] += 360
        model_data_combined = model_data_combined.assign_coords(longitude=model_longitudes)

    # Convert model time to UNIX time and ensure it's double
    model_time_seconds = np.array((model_data_combined['time'].values - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')).astype('double')
    max_model_time = model_time_seconds.max()
    min_model_time = model_time_seconds.min()  

    # Defining model_times correctly for use in interpolation
    model_times = model_time_seconds

    # Load satellite data
    satellite_data = xr.open_dataset(satellite_file)
    satellite_longitudes = satellite_data['longitude'].values
    satellite_longitudes[satellite_longitudes < 0] += 360
    satellite_times_unix = np.array((satellite_data['time'].values - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')).astype('double')
    # Apply both min and max time filters for consistency with GRIB2 section
    time_mask = (satellite_times_unix >= min_model_time) & (satellite_times_unix <= max_model_time)
    filtered_satellite_times_unix = satellite_times_unix[time_mask]
    filtered_satellite_latitudes = satellite_data['latitude'].values[time_mask]
    filtered_satellite_longitudes = satellite_longitudes[time_mask]


    # Attempt to read the satellite name from the attributes or a variable
    if 'sat_name' in satellite_data.variables:
        satellite_name = satellite_data['sat_name'].values[0]  # Assuming it's stored as a variable
    elif 'satellite_name' in satellite_data.attrs:
        satellite_name = satellite_data.attrs['satellite_name']  # Assuming it's stored as an attribute
    else:
        satellite_name = "unknown_satellite"  # Default name if not found


    # Interpolate for each variable
    interpolated_htsgw_surface_values, interpolated_htsgw_surface_times = prepare_and_interpolate(
        model_data_combined, 'HTSGW_surface', filtered_satellite_latitudes, filtered_satellite_longitudes, filtered_satellite_times_unix, model_time_seconds)

    interpolated_wind_surface_values, interpolated_wind_surface_times = prepare_and_interpolate(
        model_data_combined, 'WIND_surface', filtered_satellite_latitudes, filtered_satellite_longitudes, filtered_satellite_times_unix, model_time_seconds)

    fcst_hr2 = np.full_like(interpolated_wind_surface_values, np.nan, dtype='float64')
    valid_indices = ~np.isnan(interpolated_wind_surface_values)
    fcst_hr2[valid_indices] = (filtered_satellite_times_unix[valid_indices] - model_reference_time_unix) / 3600

    fcst_hr = np.full_like(filtered_satellite_times_unix, np.nan, dtype='float64')
    for i, value in enumerate(interpolated_wind_surface_values):
        if not np.isnan(value):
            fcst_hr[i] = filtered_satellite_times_unix[i]

    reference_time, reference_date = extract_reference_times_netcdf(model_data_files[0])

    # Extract hs and wsp from satellite data, filtered by time
    hs = satellite_data['hs'].values[time_mask]
    wsp_cal = satellite_data['wsp_cal'].values[time_mask]
    hs_cal = satellite_data['hs_cal'].values[time_mask]
    wsp = satellite_data['wsp'].values[time_mask]


    # Convert satellite time DataArray to UNIX timestamps
    satellite_time_unix = xr.DataArray(filtered_satellite_times_unix, dims=['time'], name='time')

    time_dataarray = xr.DataArray(filtered_satellite_times_unix, dims=['time'], name='time', attrs={
        'standard_name': 'time',
        'units': 'seconds since 1970-01-01 00:00:00',
        'calendar': 'standard',
        'axis': 'T'
    })

    # Create DataArrays for all data with consistent time representation and units
    hs_dataarray = xr.DataArray(hs, coords={'time': satellite_time_unix}, dims=['time'], name='hs_time_averaged_satellite')
    hs_dataarray.attrs['units'] = 'm'
    wsp_dataarray = xr.DataArray(wsp, coords={'time': satellite_time_unix}, dims=['time'], name='wsp_time_averaged_satellite')
    wsp_dataarray.attrs['units'] = 'm/s'
    interpolated_htsgw_surface = xr.DataArray(interpolated_htsgw_surface_values, coords={'time': satellite_time_unix}, dims=['time'], name='HTSGW_surface_interpolated')
    interpolated_htsgw_surface.attrs['units'] = 'm'
    interpolated_wind_surface = xr.DataArray(interpolated_wind_surface_values, coords={'time': satellite_time_unix}, dims=['time'], name='WIND_surface_interpolated')
    interpolated_wind_surface.attrs['units'] = 'm/s'
    longitude_dataarray = xr.DataArray(filtered_satellite_longitudes, coords={'time': satellite_time_unix}, dims=['time'], name='longitude')
    longitude_dataarray.attrs['units'] = 'degrees_east'
    latitude_dataarray = xr.DataArray(filtered_satellite_latitudes, coords={'time': satellite_time_unix}, dims=['time'], name='latitude')
    latitude_dataarray.attrs['units'] = 'degrees_north'
    hs_cal_dataarray = xr.DataArray(hs_cal, coords={'time': satellite_time_unix}, dims=['time'], name='hs_cal_time_averaged_satellite')
    hs_cal_dataarray.attrs['units'] = 'm'
    wsp_cal_dataarray = xr.DataArray(wsp_cal, coords={'time': satellite_time_unix}, dims=['time'], name='wsp_cal_time_averaged_satellite')
    wsp_dataarray.attrs['units'] = 'm/s'
    fcst_hr_dataarray = xr.DataArray(fcst_hr2, coords={'time': satellite_time_unix}, dims=['time'], name='fcst_time').assign_attrs(units='hours') #fcst_hr


    # Combine into a single dataset / 'time': satellite_time_unix
    interpolated_dataset = xr.Dataset({
        'time': time_dataarray,
        'model_hs': interpolated_htsgw_surface,
        'model_wnd': interpolated_wind_surface,
        'obs_hs': hs_dataarray,
        'obs_wnd': wsp_dataarray,
        'longitude': longitude_dataarray,
        'latitude': latitude_dataarray,
        'obs_hs_cal':hs_cal_dataarray,
        'obs_wnd_cal':wsp_dataarray,
        'fcst_hr': fcst_hr_dataarray
    })

    # Before saving your dataset, add the reference times as attributes
    interpolated_dataset.attrs['satellite_name'] = satellite_name
    interpolated_dataset.attrs['initial_condition_time'] = reference_time
    interpolated_dataset.attrs['initial_condition_time_units'] = "seconds since 1970-01-01 00:00:00"  # Specify units separately
    base_output_filename = "ProcessedData"
    interpolated_dataset.attrs['model_name'] = model_name

    # Save the combined dataset to a NetCDF file
    output_filename = f"{base_output_filename}_{satellite_name}_{reference_date}.nc"
    interpolated_dataset.to_netcdf(output_file, format='NETCDF4')


def spatial_mask(satellite_latitudes, satellite_longitudes,model_data):
    min_lat_val = model_data['latitude'].min().values.item()
    max_lat_val = model_data['latitude'].max().values.item()
    min_lon_val = model_data['longitude'].min().values.item()
    max_lon_val = model_data['longitude'].max().values.item()

    within_lat_bounds = (satellite_latitudes >= min_lat_val) & (satellite_latitudes <= max_lat_val)
    within_lon_bounds = (satellite_longitudes >= min_lon_val) & (satellite_longitudes <= max_lon_val)
    return within_lat_bounds & within_lon_bounds

def prepare_and_interpolate(model_data, variable_name, satellite_latitudes, satellite_longitudes, satellite_times_unix, valid_times_unix):
    mask = spatial_mask(satellite_latitudes, satellite_longitudes, model_data)
    valid_satellite_points = np.array([satellite_latitudes[mask], satellite_longitudes[mask], satellite_times_unix[mask]]).T

    variable_data = model_data[variable_name]
    # Ensure the data is in the correct dimension order
    if variable_data.dims != ('latitude', 'longitude', 'time'):
        variable_data = variable_data.transpose('latitude', 'longitude', 'time')

    model_grid = (variable_data.coords['latitude'].values,
                  variable_data.coords['longitude'].values,
                  valid_times_unix)

    interpolator = RegularGridInterpolator(model_grid, variable_data.values, bounds_error=False, fill_value=np.nan)
    interpolated_values = interpolator(valid_satellite_points)

    # Handle model times for interpolated values based on satellite times
    model_times_for_interpolated_values = np.searchsorted(valid_times_unix, satellite_times_unix[mask], side='right') - 1
    model_times_for_interpolated_values = valid_times_unix[model_times_for_interpolated_values]

    # Apply NaN where interpolated values are NaN
    nan_indices = np.isnan(interpolated_values)
    model_times_for_interpolated_values[nan_indices] = np.nan

    full_interpolated_values = np.full(len(satellite_times_unix), np.nan)
    full_interpolated_values[mask] = interpolated_values

    full_model_times = np.full(len(satellite_times_unix), np.nan)
    full_model_times[mask] = model_times_for_interpolated_values

    return full_interpolated_values, full_model_times



def main():


    ap = argparse.ArgumentParser()
    ap.add_argument('-t', '--typefile', help="Type of Model File, 'grib2' or 'nc'",required=True)
    ap.add_argument('-d', '--datadir', help="Data Directory for Model Files", required=True)
    ap.add_argument('-p', '--pattern', help="Pattern of Model Files", required=True)
    ap.add_argument('-s', '--satfile', help="Satellite File", required=True)
    ap.add_argument('-o', '--outdir', help="Directory Path for Output", required=True)
    ap.add_argument('-f', '--fileout', help="Name of Output File", required=True)
    ap.add_argument('-m', '--model', help="String Identifier of Model", required=True)

    MyArgs = ap.parse_args()

    file_type = MyArgs.typefile
    data_directory = MyArgs.datadir
    data_pattern = MyArgs.pattern
    satellite_file = MyArgs.satfile
    path_out = MyArgs.outdir
    outfilename = MyArgs.fileout
    model_name = MyArgs.model

    print(f"File type: {file_type}")
    print(f"Data directory: {data_directory}")
    print(f"Data pattern: {data_pattern}")
    print(f"Satellite file: {satellite_file}")
    print(f"Output directory: {path_out}")
    print(f"Output file: {outfilename}")
    print(f"Model name: {model_name}")

    #create path_out directory if it does not exist:
    if not os.path.isdir(path_out):
        os.makedirs(path_out)

    output_file = path_out + "/" + outfilename

    # Call the appropriate function based on file type
    if file_type == 'grib2':
        interpolate_grib2(data_directory, data_pattern, satellite_file, output_file, model_name)
    elif file_type == 'nc':
        interpolate_netcdf(data_directory,data_pattern, satellite_file, output_file, model_name)
    else:
        raise ValueError("Invalid file type specified. Please use 'grib2' or 'nc'.")

if __name__ == "__main__":
    main()

