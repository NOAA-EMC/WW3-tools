"""
This script processes oceanographic model data and satellite time avaraged observations, aligning and interpolating them for comparison and analysis. It is specifically designed to work with gridded model output and satellite data, focusing on significant wave height (HTSGW,swh) and wind speed (WIND,ws) parameters.

Key Steps:
1. Data Loading: The script loads model data from a specified directory, where the files are expected to follow a naming pattern. It also loads satellite data from a given file.

2. Time Alignment: Both model and satellite data times are converted to UNIX time for consistent time referencing. The script extracts forecast hours from model filenames and aligns them with satellite observation times.

3. Data Preparation and Interpolation: The model data is prepared by applying an ocean mask to handle land values. Interpolation is then performed using a regular grid interpolator, aligning the model data with the satellite data points.

4. Data Merging: The interpolated model data and the satellite observations are combined into a single xarray Dataset. This includes significant wave height and wind speed parameters from both sources, as well as the aligned forecast hours.

5. Output: The final combined dataset is saved as a NetCDF file, which can be used for further analysis, visualization, or comparison studies.

Requirements:
- xarray: For handling NetCDF files and data arrays.
- numpy: For numerical operations.
- glob: For file path pattern matching.
- re: For regular expression operations.
- scipy.interpolate: For data interpolation on regular grids.
- pandas: For data manipulation and analysis.
-cfgrib: cfgrib is a Python library that provides an interface for reading GRIB (GRIdded Binary) files, which are a common format for storing and distributing
meteorological and oceanographic model data.

The script is intended for oceanographic data analysts and researchers who need to compare model output with satellite measurements for validation or study purposes.
In this python script:
file_type = sys.argv[1].lower()  # Expecting 'grib2' or 'nc'
data_directory = sys.argv[2]
data_pattern = sys.argv[3]
satellite_file = sys.argv[4]
output_file = sys.argv[5]
model_name = sys.argv[6]  # New argument for satellite name
ic_time = sys.argv[7]  # New argument for IC time

are pathes for the files thay one can define in the job script.
model_data_directory: The directory where model data files are located. For example, /scratch2/NCEPDEV/marine/Jessica.Meixner/Data/HR1/Hurricane/gfs.20200919/00/wave/gridded/.
model_data_pattern: A pattern to match specific model data files within the directory. Example: 'gfswave.t00z.global.0p25.f*.grib2.nc' which selects files with a certain naming convention.
satellite_file: The path to the satellite data file. Example: './AltimeterAlongTrack_ww3tools_JASON3_2020091901to2020092922.nc'.
output_file: The path for the output file where processed results will be saved. Example: './WW3-Altimeter_interpolated_20200919.nc'.


# Function Definitions:

def interpolate_grib2(grib_data_directory, grib_data_pattern, satellite_file, output_file, satellite_name, ic_time):

Handles the interpolation of GRIB2 format model data.

This function searches for GRIB2 files in a specified directory, aligns them with satellite observation times,
and interpolates the model data to the satellite data points. The output is a combined dataset in NetCDF format.

Parameters:
- grib_data_directory (str): Path to the directory containing the GRIB2 files.
- grib_data_pattern (str): Pattern to match specific GRIB2 files within the directory.
- satellite_file (str): Path to the satellite data file.
- output_file (str): Path for the output NetCDF file.

The function processes each GRIB file, extracting significant wave height and wind speed data, aligning them in time
with satellite observations, and performing spatial interpolation. The result is saved as a NetCDF file.
# Function code...

def interpolate_netcdf(model_data_directory, model_data_pattern, satellite_file, output_file):

Handles the interpolation of NetCDF format model data.

Similar to 'interpolate_grib2', this function works with NetCDF format model data. It loads the data, aligns it
with satellite observations, and performs the necessary interpolations. The result is a merged dataset in NetCDF format.

Parameters:
- model_data_directory (str): Path to the directory containing the NetCDF files.
- model_data_pattern (str): Pattern to match specific NetCDF files within the directory.
- satellite_file (str): Path to the satellite data file.
- output_file (str): Path for the output NetCDF file.

This function reads the NetCDF files, applies a spatial mask to handle ocean and land values, and interpolates the data
to align with satellite data points. The interpolated dataset includes parameters like significant wave height and wind speed.


Author: Ghazal Mohammadpour
email: ghazal.mohammadpour@noaa.gov

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

#GRIB2 interpolator

# Function to map forecast hours to satellite times
def map_fcst_hours_to_satellite(fcst_hours, model_times_unix, satellite_times_unix):
    df = pd.DataFrame({'fcst_hr': fcst_hours}, index=model_times_unix)
    df = df.reindex(satellite_times_unix, method='nearest')
    return df['fcst_hr'].values


def interpolate_grib2(data_directory, data_pattern, satellite_file, output_file, model_name, ic_time):
    full_path = os.path.join(data_directory, data_pattern)
    print("Full path for GRIB files:", full_path)

    found_files = glob.glob(full_path)
    if not found_files:
        print("No files found at:", full_path)
        raise ValueError("No GRIB files found in the specified directory matching the pattern.")

    found_files.sort()
    fcst_hours = [int(re.search(r'f(\d+)', file).group(1)) for file in found_files]
    #fcst_hours_filtered = [hour for hour in fcst_hours if 0 <= hour <= 24]
    fcst_hours_filtered = fcst_hours

    # Extract satellite name from the satellite file name
#    satellite_name_search = re.search(r'_([^-]+)_', satellite_file)
#    satellite_name = satellite_name_search.group(1) if satellite_name_search else "unknown_satellite"

    satellite_name_search = re.search(r'ww3tools_([^_]+)', satellite_file)
    satellite_name = satellite_name_search.group(1) if satellite_name_search else "unknown_satellite"


    filtered_grib_files = [found_files[i] for i, hour in enumerate(fcst_hours) if hour in fcst_hours_filtered]
    model_data_list, valid_times_unix = process_grib_files(filtered_grib_files)
    model_data_combined = xr.concat(model_data_list, dim='time')

    satellite_data = xr.open_dataset(satellite_file)
    satellite_times_unix = np.array((satellite_data['time'].values - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')).astype('double')

    min_model_time = min(valid_times_unix)
    max_model_time = max(valid_times_unix)

    time_mask = (satellite_times_unix >= min_model_time) & (satellite_times_unix <= max_model_time)
    filtered_satellite_times_unix = satellite_times_unix[time_mask]

    mapped_fcst_hours = map_fcst_hours_to_satellite(fcst_hours, valid_times_unix, filtered_satellite_times_unix)


    # Extract additional satellite variables
    hs = satellite_data['hs'].values[time_mask]
    wsp_cal = satellite_data['wsp_cal'].values[time_mask]
    satellite_longitudes = satellite_data['longitude'].values[time_mask]
    satellite_latitudes = satellite_data['latitude'].values[time_mask]
    hs_cal = satellite_data['hs_cal'].values[time_mask]
    wsp = satellite_data['wsp'].values[time_mask]

   # Modify here to add time attributes
    time_dataarray = xr.DataArray(filtered_satellite_times_unix, dims=['time'], name='time', attrs={
        'standard_name': 'time',
        'units': 'seconds since 1970-01-01 00:00:00',
        'calendar': 'standard',
        'axis': 'T'
    })

    interpolated_ws_values = prepare_and_interpolate(model_data_combined, 'wind_speed', satellite_latitudes, satellite_longitudes, filtered_satellite_times_unix, valid_times_unix)
    interpolated_swh_values = prepare_and_interpolate(model_data_combined, 'significant_wave_height', satellite_latitudes, satellite_longitudes, filtered_satellite_times_unix, valid_times_unix)

    # Create DataArrays for all variables
    interpolated_ws_dataarray = xr.DataArray(interpolated_ws_values, coords={'time': filtered_satellite_times_unix}, dims=['time'], name='ws_interpolated')
    interpolated_swh_dataarray = xr.DataArray(interpolated_swh_values, coords={'time': filtered_satellite_times_unix}, dims=['time'], name='swh_interpolated')
    hs_dataarray = xr.DataArray(hs, coords={'time': filtered_satellite_times_unix}, dims=['time'], name='hs')
    wsp_cal_dataarray = xr.DataArray(wsp_cal, coords={'time': filtered_satellite_times_unix}, dims=['time'], name='wsp_cal')
    longitude_dataarray = xr.DataArray(satellite_longitudes, coords={'time': filtered_satellite_times_unix}, dims=['time'], name='longitude')
    latitude_dataarray = xr.DataArray(satellite_latitudes, coords={'time': filtered_satellite_times_unix}, dims=['time'], name='latitude')
    hs_cal_dataarray = xr.DataArray(hs_cal, coords={'time': filtered_satellite_times_unix}, dims=['time'], name='hs_cal')
    wsp_dataarray = xr.DataArray(wsp, coords={'time': filtered_satellite_times_unix}, dims=['time'], name='wsp')
    fcst_hours_dataarray = xr.DataArray(mapped_fcst_hours, coords={'time': filtered_satellite_times_unix}, dims=['time'], name='fcst_hr')
    fcst_hours_dataarray.attrs['units'] = 'hours'

    # Combine all DataArrays into a single Dataset
    interpolated_dataset = xr.Dataset({
        'time':time_dataarray,
        'ws_interpolated': interpolated_ws_dataarray,
        'swh_interpolated': interpolated_swh_dataarray,
        'hs': hs_dataarray,
        'wsp_cal': wsp_cal_dataarray,
        'longitude': longitude_dataarray,
        'latitude': latitude_dataarray,
        'hs_cal':hs_cal_dataarray,
        'wsp':wsp_dataarray,
        'fcst_hr': fcst_hours_dataarray
    })

    interpolated_dataset.attrs['satellite_name']= satellite_name
    interpolated_dataset.attrs['IC'] = ic_time
    interpolated_dataset.attrs['model_name'] = model_name
    interpolated_dataset.to_netcdf(output_file, format='NETCDF4')

def process_grib_files(files):
    model_data_list = []
    valid_times_unix = []

    for file in files:
        ds = cfgrib.open_datasets(file, backend_kwargs={'indexpath': ''})
        if len(ds) < 2:
            raise ValueError(f"Required datasets not found in the file {file}")

        ws = ds[1]['ws'] if 'ws' in ds[1] else xr.full_like(ds[1]['latitude'], np.nan, dtype=np.float32)
        swh = ds[1]['swh'] if 'swh' in ds[1] else xr.full_like(ds[1]['latitude'], np.nan, dtype=np.float32)
        longitude = ds[1]['longitude']
        latitude = ds[1]['latitude']
        time = ds[1]['time']
        step = ds[1]['step']

        valid_time = ds[1].valid_time.values
        valid_time_unix = (valid_time - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')
        valid_times_unix.append(valid_time_unix)

        data = xr.Dataset({
            'wind_speed': ws,
            'significant_wave_height': swh,
            'longitude': longitude,
            'latitude': latitude,
            'time': time,
            'step': step
        })

        model_data_list.append(data)

    return model_data_list, valid_times_unix

def prepare_and_interpolate(model_data, variable_name, satellite_latitudes, satellite_longitudes, satellite_times, valid_times_unix):
    mask = spatial_mask(satellite_latitudes, satellite_longitudes, model_data)
    valid_satellite_points = np.array([satellite_latitudes[mask], satellite_longitudes[mask], satellite_times[mask]]).T

    variable_data = model_data[variable_name]
    if variable_data.dims != ('latitude', 'longitude', 'time'):
        variable_data = variable_data.transpose('latitude', 'longitude', 'time')

    ocean_mask = ~variable_data.isnull()
    variable_data_masked = variable_data.where(ocean_mask, np.nan)

    model_grid = (variable_data.coords['latitude'].values,
                  variable_data.coords['longitude'].values,
                  valid_times_unix)

    if valid_satellite_points.size == 0:
        return np.full((len(satellite_times),), np.nan)

    interpolator = RegularGridInterpolator(model_grid, variable_data_masked.values, bounds_error=False, fill_value=np.nan)
    interpolated_values = interpolator(valid_satellite_points)

    full_interpolated_values = np.full((len(satellite_times),), np.nan)
    full_interpolated_values[mask] = interpolated_values

    return full_interpolated_values

def spatial_mask(satellite_latitudes, satellite_longitudes, model_data):
    min_latitude = model_data.latitude.min().values.item()
    max_latitude = model_data.latitude.max().values.item()
    min_longitude = model_data.longitude.min().values.item()
    max_longitude = model_data.longitude.max().values.item()

    within_lat_bounds = (satellite_latitudes >= min_latitude) & (satellite_latitudes <= max_latitude)
    within_lon_bounds = (satellite_longitudes >= min_longitude) & (satellite_longitudes <= max_longitude)
    return within_lat_bounds & within_lon_bounds



#NETCDF interpolator


def interpolate_netcdf(model_data_directory, model_data_pattern, satellite_file, output_file, model_name, ic_time):
    # Using glob to find all files in the directory that match the pattern
    model_data_files = glob.glob(model_data_directory + model_data_pattern)

    # Check if any files are found
    if not model_data_files:
        raise ValueError("No files found in the specified directory matching the pattern.")

    # Sort the files to ensure they are in the correct order
    model_data_files.sort()

    # Extract forecast hours from file names and convert to integers
    fcst_hours = [int(re.search(r'f(\d+)', file).group(1)) for file in model_data_files]
    #fcst_hours_filtered = [hour for hour in fcst_hours if 144 <= hour <= 168]
    fcst_hours_filtered = fcst_hours

    # Extract satellite name from the satellite file name
    satellite_name_search = re.search(r'ww3tools_([^_]+)', satellite_file)
    satellite_name = satellite_name_search.group(1) if satellite_name_search else "unknown_satellite"

    # Load and concatenate the model data files
    model_data_list = [xr.open_dataset(file) for file in model_data_files if int(re.search(r'f(\d+)', file).group(1)) in fcst_hours_filtered]
    model_data_combined = xr.concat(model_data_list, dim='time')

    # Convert longitudes from -180 to 180 to 0 to 360, if necessary
    model_longitudes = model_data_combined['longitude'].values
    if model_longitudes.min() < 0:
        model_longitudes[model_longitudes < 0] += 360
        model_data_combined = model_data_combined.assign_coords(longitude=model_longitudes)

    # Convert model time to UNIX time and ensure it's double
    model_time_seconds = np.array((model_data_combined['time'].values - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')).astype('double')
    max_model_time = model_time_seconds.max()

    # Load satellite data
    satellite_data = xr.open_dataset(satellite_file)
    satellite_longitudes = satellite_data['longitude'].values
    satellite_longitudes[satellite_longitudes < 0] += 360
    satellite_times_unix = np.array((satellite_data['time'].values - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')).astype('double')
    filtered_satellite_times_unix = satellite_times_unix[satellite_times_unix <= max_model_time]
    filtered_satellite_latitudes = satellite_data['latitude'].values[:len(filtered_satellite_times_unix)]
    filtered_satellite_longitudes = satellite_longitudes[:len(filtered_satellite_times_unix)]

    # Map forecast hours to satellite times
    mapped_fcst_hours = map_fcst_hours_to_satellite(fcst_hours, model_time_seconds, satellite_times_unix)

    # Define the Model's Spatial Bounds
    min_latitude = model_data_combined['latitude'].min()
    max_latitude = model_data_combined['latitude'].max()
    min_longitude = model_data_combined['longitude'].min()
    max_longitude = model_data_combined['longitude'].max()

    # Function to create and apply a spatial mask
    def spatial_mask(satellite_latitudes, satellite_longitudes):
        min_lat_val = min_latitude.values.item()
        max_lat_val = max_latitude.values.item()
        min_lon_val = min_longitude.values.item()
        max_lon_val = max_longitude.values.item()

        within_lat_bounds = (satellite_latitudes >= min_lat_val) & (satellite_latitudes <= max_lat_val)
        within_lon_bounds = (satellite_longitudes >= min_lon_val) & (satellite_longitudes <= max_lon_val)
        return within_lat_bounds & within_lon_bounds

    # Function to prepare and interpolate model data
    def prepare_and_interpolate(model_data, variable_name, satellite_latitudes, satellite_longitudes, satellite_times_unix):
        mask = spatial_mask(satellite_latitudes, satellite_longitudes)
        valid_satellite_points = np.array([satellite_latitudes[mask], satellite_longitudes[mask], satellite_times_unix[mask]]).T

        variable_data = model_data[variable_name].transpose('latitude', 'longitude', 'time')
        ocean_mask = ~variable_data.isnull()
        variable_data_masked = variable_data.where(ocean_mask, np.nan)
        model_grid = (variable_data.coords['latitude'].values, variable_data.coords['longitude'].values, model_time_seconds)
        interpolator = RegularGridInterpolator(model_grid, variable_data_masked.values, bounds_error=False, fill_value=np.nan)

        interpolated_values = interpolator(valid_satellite_points)
        full_interpolated_values = np.full(satellite_latitudes.shape, np.nan)
        full_interpolated_values[mask] = interpolated_values

        return full_interpolated_values

    # Perform interpolation using the modified function
    interpolated_htsgw_surface_values = prepare_and_interpolate(model_data_combined, 'HTSGW_surface', filtered_satellite_latitudes, filtered_satellite_longitudes, filtered_satellite_times_unix)
    interpolated_wind_surface_values = prepare_and_interpolate(model_data_combined, 'WIND_surface', filtered_satellite_latitudes, filtered_satellite_longitudes, filtered_satellite_times_unix)

    # Extract hs and wsp from satellite data, filtered by time
    hs = satellite_data['hs'].sel(time=satellite_data['time'].values[satellite_times_unix <= max_model_time])
    wsp_cal = satellite_data['wsp_cal'].sel(time=satellite_data['time'].values[satellite_times_unix <= max_model_time])  #I dont added _cal to the name
    hs_cal = satellite_data['hs_cal'].sel(time=satellite_data['time'].values[satellite_times_unix <= max_model_time])
    wsp = satellite_data['wsp'].sel(time=satellite_data['time'].values[satellite_times_unix <= max_model_time])

    # Convert satellite time DataArray to UNIX timestamps
    satellite_time_unix = xr.DataArray(filtered_satellite_times_unix, dims=['time'], name='time')

    time_dataarray = xr.DataArray(filtered_satellite_times_unix, dims=['time'], name='time', attrs={
        'standard_name': 'time',
        'units': 'seconds since 1970-01-01 00:00:00',
        'calendar': 'standard',
        'axis': 'T'
    })

    # Create DataArrays for all data with consistent time representation and units
    hs_dataarray = xr.DataArray(hs.values, coords={'time': satellite_time_unix}, dims=['time'], name='hs_time_averaged_satellite')
    hs_dataarray.attrs['units'] = 'm'
    wsp_dataarray = xr.DataArray(wsp.values, coords={'time': satellite_time_unix}, dims=['time'], name='wsp_time_averaged_satellite')
    wsp_dataarray.attrs['units'] = 'm/s'
    interpolated_htsgw_surface = xr.DataArray(interpolated_htsgw_surface_values, coords={'time': satellite_time_unix}, dims=['time'], name='HTSGW_surface_interpolated')
    interpolated_htsgw_surface.attrs['units'] = 'm'
    interpolated_wind_surface = xr.DataArray(interpolated_wind_surface_values, coords={'time': satellite_time_unix}, dims=['time'], name='WIND_surface_interpolated')
    interpolated_wind_surface.attrs['units'] = 'm/s'
    longitude_dataarray = xr.DataArray(filtered_satellite_longitudes, coords={'time': satellite_time_unix}, dims=['time'], name='longitude')
    longitude_dataarray.attrs['units'] = 'degrees_east'
    latitude_dataarray = xr.DataArray(filtered_satellite_latitudes, coords={'time': satellite_time_unix}, dims=['time'], name='latitude')
    latitude_dataarray.attrs['units'] = 'degrees_north'
    hs_cal_dataarray = xr.DataArray(hs_cal.values, coords={'time': satellite_time_unix}, dims=['time'], name='hs_cal_time_averaged_satellite')
    hs_cal_dataarray.attrs['units'] = 'm'
    wsp_cal_dataarray = xr.DataArray(wsp_cal.values, coords={'time': satellite_time_unix}, dims=['time'], name='wsp_cal_time_averaged_satellite')
    wsp_dataarray.attrs['units'] = 'm/s'
    fcst_hours_dataarray = xr.DataArray(mapped_fcst_hours, coords={'time': satellite_times_unix}, dims=['time'], name='fcst_hr')
    fcst_hours_dataarray.attrs['units'] = 'hours'

    # Combine into a single dataset / 'time': satellite_time_unix
    interpolated_dataset = xr.Dataset({
        'time': time_dataarray,
        'HTSGW_surface_interpolated': interpolated_htsgw_surface,
        'WIND_surface_interpolated': interpolated_wind_surface,
        'hs_time_averaged_satellite': hs_dataarray,
        'wsp_time_averaged_satellite': wsp_dataarray,
        'longitude': longitude_dataarray,
        'latitude': latitude_dataarray,
        'hs_cal_time_averaged_satellite':hs_cal_dataarray,
        'wsp_cal_time_averaged_satellite':wsp_dataarray,
        'fcst_hr': fcst_hours_dataarray
    })


    interpolated_dataset.attrs['satellite_name'] = satellite_name
    interpolated_dataset.attrs['IC'] = ic_time
    interpolated_dataset.attrs['model_name'] = model_name
    # Save the combined dataset to a NetCDF file
    interpolated_dataset.to_netcdf(output_file, format='NETCDF4')








def main():
    # Command-line arguments
    file_type = sys.argv[1].lower()  # Expecting 'grib2' or 'nc'
    data_directory = sys.argv[2]
    data_pattern = sys.argv[3]
    satellite_file = sys.argv[4]
    output_file = sys.argv[5]
    model_name = sys.argv[6]  # New argument for satellite name
    ic_time = sys.argv[7]  # New argument for IC time

    print(f"File type: {file_type}")
    print(f"Data directory: {data_directory}")
    print(f"Data pattern: {data_pattern}")
    print(f"Satellite file: {satellite_file}")
    print(f"Output file: {output_file}")
    print(f"model name: {model_name}")
    print(f"IC: {ic_time}")

    # Call the appropriate function based on file type
    if file_type == 'grib2':
        interpolate_grib2(data_directory, data_pattern, satellite_file, output_file, model_name, ic_time)
    elif file_type == 'nc':
        interpolate_netcdf(data_directory, data_pattern, satellite_file, output_file, model_name, ic_time)
    else:
        raise ValueError("Invalid file type specified. Please use 'grib2' or 'nc'.")

if __name__ == "__main__":
    main()



