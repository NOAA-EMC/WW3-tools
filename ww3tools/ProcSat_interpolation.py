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

are pathes for the files thay one can define in the job script.
model_data_directory: The directory where model data files are located. For example, /scratch2/NCEPDEV/marine/Jessica.Meixner/Data/HR1/Hurricane/gfs.20200919/00/wave/gridded/.
model_data_pattern: A pattern to match specific model data files within the directory. Example: 'gfswave.t00z.global.0p25.f*.grib2.nc' which selects files with a certain naming convention.
satellite_file: The path to the satellite data file. Example: './AltimeterAlongTrack_ww3tools_JASON3_2020091901to2020092922.nc'.
output_file: The path for the output file where processed results will be saved. Example: './WW3-Altimeter_interpolated_20200919.nc'.


# Function Definitions:

def interpolate_grib2(grib_data_directory, grib_data_pattern, satellite_file, output_file):
    """
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
    """
    # Function code...

def interpolate_netcdf(model_data_directory, model_data_pattern, satellite_file, output_file):
    """
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
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
import cfgrib

def interpolate_grib2(grib_data_directory, grib_data_pattern, satellite_file, output_file):
    # Using glob to find all GRIB files in the directory that match the pattern
    grib_data_files = glob.glob(grib_data_directory + grib_data_pattern)

    # Check if any files are found
    if not grib_data_files:
        raise ValueError("No GRIB files found in the specified directory matching the pattern.")

    # Sort the files and extract forecast hours
    grib_data_files.sort()
    fcst_hours = [int(re.search(r'f(\d+)', file).group(1)) for file in grib_data_files]

    # Filter forecast hours to include only 0 to 24 hours
    fcst_hours_filtered = [hour for hour in fcst_hours if 0 <= hour <= 24]

    # Function to read and process each GRIB file
    def process_grib_files(files):
        model_data_list = []
        valid_times_unix = []

        for file in files:
            ds = cfgrib.open_datasets(file)
            if len(ds) < 2:
                raise ValueError(f"Required datasets not found in the file {file}")

            # Extract variables from the second dataset
            ws = ds[1]['ws']  # example variable name
            swh = ds[1]['swh']  # example variable name
            longitude = ds[1]['longitude']
            latitude = ds[1]['latitude']
            time = ds[1]['time']
            step = ds[1]['step']

            # Extract valid_time and convert to UNIX time
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

    # Process all GRIB files with filtered forecast hours
    filtered_grib_files = [grib_data_files[i] for i, hour in enumerate(fcst_hours) if hour in fcst_hours_filtered]
    model_data_list, valid_times_unix = process_grib_files(filtered_grib_files)

    # Concatenate the list of datasets into a single dataset
    model_data_combined = xr.concat(model_data_list, dim='time')

    # Load satellite data
    satellite_data = xr.open_dataset(satellite_file)

    # Convert satellite time to UNIX time
    satellite_times_unix = np.array((satellite_data['time'].values - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')).astype('double')

    # Define time limits based on the model's valid times
    min_model_time = min(valid_times_unix)
    max_model_time = max(valid_times_unix)

    # Filter satellite times to only include those within the range of model times
    time_mask = (satellite_times_unix >= min_model_time) & (satellite_times_unix <= max_model_time)
    filtered_satellite_times_unix = satellite_times_unix[time_mask]
    filtered_satellite_latitudes = satellite_data['latitude'].values[time_mask]
    filtered_satellite_longitudes = satellite_data['longitude'].values[time_mask]

    # Spatial masking function
    def spatial_mask(satellite_latitudes, satellite_longitudes):
        min_latitude = model_data_combined.latitude.min().values.item()
        max_latitude = model_data_combined.latitude.max().values.item()
        min_longitude = model_data_combined.longitude.min().values.item()
        max_longitude = model_data_combined.longitude.max().values.item()

        within_lat_bounds = (satellite_latitudes >= min_latitude) & (satellite_latitudes <= max_latitude)
        within_lon_bounds = (satellite_longitudes >= min_longitude) & (satellite_longitudes <= max_longitude)
        return within_lat_bounds & within_lon_bounds

    # Modified interpolation function
    def prepare_and_interpolate(model_data, variable_name, satellite_latitudes, satellite_longitudes, satellite_times):
        mask = spatial_mask(satellite_latitudes, satellite_longitudes)
        valid_satellite_points = np.array([satellite_latitudes[mask], satellite_longitudes[mask], satellite_times[mask]]).T

        variable_data = model_data[variable_name]
        if variable_data.dims != ('latitude', 'longitude', 'time'):
            variable_data = variable_data.transpose('latitude', 'longitude', 'time')

        ocean_mask = ~variable_data.isnull()
        variable_data_masked = variable_data.where(ocean_mask, np.nan)

        model_grid = (variable_data.coords['latitude'].values,
                      variable_data.coords['longitude'].values,
                      valid_times_unix)

        interpolator = RegularGridInterpolator(model_grid, variable_data_masked.values, bounds_error=False, fill_value=np.nan)
        interpolated_values = interpolator(valid_satellite_points)

        return interpolated_values

    # Perform interpolation for each variable
    interpolated_ws_values = prepare_and_interpolate(model_data_combined, 'wind_speed', filtered_satellite_latitudes, filtered_satellite_longitudes, filtered_satellite_times_unix)
    interpolated_swh_values = prepare_and_interpolate(model_data_combined, 'significant_wave_height', filtered_satellite_latitudes, filtered_satellite_longitudes, filtered_satellite_times_unix)

    # Create DataArrays for interpolated values
    interpolated_ws_dataarray = xr.DataArray(interpolated_ws_values, coords={'time': filtered_satellite_times_unix}, dims=['time'], name='ws_interpolated')
    interpolated_swh_dataarray = xr.DataArray(interpolated_swh_values, coords={'time': filtered_satellite_times_unix}, dims=['time'], name='swh_interpolated')

    # Combine into a single dataset
    interpolated_dataset = xr.Dataset({
        'ws_interpolated': interpolated_ws_dataarray,
        'swh_interpolated': interpolated_swh_dataarray
    })

    # Adding satellite variables
    interpolated_dataset['hs'] = xr.DataArray(
        satellite_data['hs'].values[time_mask],
        coords={'time': filtered_satellite_times_unix},
        dims=['time'],
        attrs={
            '_FillValue': np.nan,
            'units': 'm'
        }
    )

    interpolated_dataset['wsp_cal'] = xr.DataArray(
        satellite_data['wsp_cal'].values[time_mask],
        coords={'time': filtered_satellite_times_unix},
        dims=['time'],
        attrs={
            '_FillValue': np.nan,
            'units': 'm/s'
        }
    )

    interpolated_dataset['longitude'] = xr.DataArray(
        filtered_satellite_longitudes,
        coords={'time': filtered_satellite_times_unix},
        dims=['time'],
        attrs={
            '_FillValue': np.nan,
            'units': 'degrees_east'
        }
    )

    interpolated_dataset['latitude'] = xr.DataArray(
        filtered_satellite_latitudes,
        coords={'time': filtered_satellite_times_unix},
        dims=['time'],
        attrs={
            '_FillValue': np.nan,
            'units': 'degrees_north'
        }
    )

    # Save the combined dataset to a NetCDF file
    interpolated_dataset.to_netcdf(output_file, format='NETCDF4')

# The rest of your script (including the main function and any other code) goes here

def interpolate_netcdf(model_data_directory, model_data_pattern, satellite_file, output_file):
    # Using glob to find all files in the directory that match the pattern
    model_data_files = glob.glob(model_data_directory + model_data_pattern)

    # Check if any files are found
    if not model_data_files:
        raise ValueError("No files found in the specified directory matching the pattern.")

    # Sort the files to ensure they are in the correct order
    model_data_files.sort()

    # Extract forecast hours from file names and convert to integers
    fcst_hours = [int(re.search(r'f(\d+)', file).group(1)) for file in model_data_files]
    fcst_hours_filtered = [hour for hour in fcst_hours if 144 <= hour <= 168]

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
    wsp = satellite_data['wsp_cal'].sel(time=satellite_data['time'].values[satellite_times_unix <= max_model_time])

    # Convert satellite time DataArray to UNIX timestamps
    satellite_time_unix = xr.DataArray(filtered_satellite_times_unix, dims=['time'], name='time')

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

    # Combine into a single dataset
    interpolated_dataset = xr.Dataset({
        'time': satellite_time_unix,
        'HTSGW_surface_interpolated': interpolated_htsgw_surface,
        'WIND_surface_interpolated': interpolated_wind_surface,
        'hs_time_averaged_satellite': hs_dataarray,
        'wsp_time_averaged_satellite': wsp_dataarray,
        'longitude': longitude_dataarray,
        'latitude': latitude_dataarray
    })

    # Save the combined dataset to a NetCDF file
    interpolated_dataset.to_netcdf(output_file, format='NETCDF4')





def main():
    # Command-line arguments
    file_type = sys.argv[1].lower()  # Expecting 'grib2' or 'nc'
    data_directory = sys.argv[2]
    data_pattern = sys.argv[3]
    satellite_file = sys.argv[4]
    output_file = sys.argv[5]

    # Call the appropriate function based on file type
    if file_type == 'grib2':
        interpolate_grib2(data_directory, data_pattern, satellite_file, output_file)
    elif file_type == 'nc':
        interpolate_netcdf(data_directory, data_pattern, satellite_file, output_file)
    else:
        raise ValueError("Invalid file type specified. Please use 'grib2' or 'nc'.")

if __name__ == "__main__":
    main()
