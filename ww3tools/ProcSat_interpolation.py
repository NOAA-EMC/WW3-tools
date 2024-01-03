"""
This script processes oceanographic model data and satellite time avaraged observations, aligning and interpolating them for comparison and analysis. It is specifically designed to work with gridded model output and satellite data, focusing on significant wave height (HTSGW) and wind speed (WIND) parameters.

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

The script is intended for oceanographic data analysts and researchers who need to compare model output with satellite measurements for validation or study purposes.
In this python script
model_data_directory = sys.argv[1]
model_data_pattern = sys.argv[2]
satellite_file = sys.argv[3]
output_file = sys.argv[4]
are pathes for the files thay one can define in the job script.
model_data_directory: The directory where model data files are located. For example, /scratch2/NCEPDEV/marine/Jessica.Meixner/Data/HR1/Hurricane/gfs.20200919/00/wave/gridded/.
model_data_pattern: A pattern to match specific model data files within the directory. Example: 'gfswave.t00z.global.0p25.f*.grib2.nc' which selects files with a certain naming convention.
satellite_file: The path to the satellite data file. Example: './AltimeterAlongTrack_ww3tools_JASON3_2020091901to2020092922.nc'.
output_file: The path for the output file where processed results will be saved. Example: './WW3-Altimeter_interpolated_20200919.nc'.


Author: Ghazal Mohammadpour
email: ghazal.mohammadpour@noaa.gov

"""

import xarray as xr
import numpy as np
import glob
import re
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
import sys

# Command-line arguments for input and output paths
model_data_directory = sys.argv[1]
model_data_pattern = sys.argv[2]
satellite_file = sys.argv[3]
output_file = sys.argv[4]

# Using glob to find all files in the directory that match the pattern
model_data_files = glob.glob(model_data_directory + model_data_pattern)

# Check if any files are found
if not model_data_files:
    raise ValueError("No files found in the specified directory matching the pattern.")

# Sort the files to ensure they are in the correct order
model_data_files.sort()

# Extract forecast hours from file names and convert to integers
fcst_hours = [int(re.search(r'f(\d+)', file).group(1)) for file in model_data_files]

# Include all forecast hours without filtering

fcst_hours_filtered = fcst_hours

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

# Determine the maximum model time
max_model_time = model_time_seconds.max()

# Load satellite data
satellite_data = xr.open_dataset(satellite_file)

# Convert satellite longitude to 0-360 format if necessary
satellite_longitudes = satellite_data['longitude'].values
satellite_longitudes[satellite_longitudes < 0] += 360

# Convert satellite time to UNIX time and ensure it's double
satellite_times_unix = np.array((satellite_data['time'].values - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')).astype('double')

# Filter satellite times to include only those within the model data time range
filtered_satellite_times_unix = satellite_times_unix[satellite_times_unix <= max_model_time]

# Filter satellite latitude and longitude arrays to match the filtered satellite times
filtered_satellite_latitudes = satellite_data['latitude'].values[:len(filtered_satellite_times_unix)]
filtered_satellite_longitudes = satellite_longitudes[:len(filtered_satellite_times_unix)]

# Function to prepare and interpolate model data
def prepare_and_interpolate(model_data, variable_name, satellite_latitudes, satellite_longitudes, satellite_times_unix):
    variable_data = model_data[variable_name].transpose('latitude', 'longitude', 'time')
    ocean_mask = ~variable_data.isnull()
    variable_data_masked = variable_data.where(ocean_mask, np.nan)
    model_grid = (variable_data.coords['latitude'].values, variable_data.coords['longitude'].values, model_time_seconds)
    interpolator = RegularGridInterpolator(model_grid, variable_data_masked.values, bounds_error=False, fill_value=np.nan)
    satellite_points = np.array([satellite_latitudes, satellite_longitudes, satellite_times_unix]).T
    return interpolator(satellite_points)

# Perform interpolation using the filtered satellite data
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

