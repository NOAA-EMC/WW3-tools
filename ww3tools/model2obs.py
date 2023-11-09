"""
This script is designed to interpolate model data onto the locations of satellite observations and combine them into a single netCDF file. Users should provide their satellite observation file and a directory containing model data files in netCDF format. The output is a netCDF file containing the satellite data along with the interpolated model data.

Steps for Users:

1. Configuration File: The user needs to create a YAML configuration file named config.yaml that includes paths to the input files and output file, as well as the interpolation method. Here's an example configuration:

input:
  obs_file: ./AltimeterAlongTrack_ww3tools_JASON3_2018010101to2018033123.nc
  model_directory: /path/to/model/data/
output:
  file: ./output_file.nc
interpolation:
  method: nearest  # Options: 'linear' or 'nearest'

Make sure to replace the paths and filenames with the actual paths to your data.

2. Input Data: Place your satellite observation .nc file and your model data directory as specified in your config.yaml.

3. Output Data: Define the desired path for the output .nc file in config.yaml. The script will generate this file with interpolated data.

4. Interpolation Method: Choose your interpolation method (nearest or linear) in config.yaml. 'Nearest' will use the nearest neighbor method, which is fast and suitable for scattered data. 'Linear' will perform linear interpolation, which is more accurate but computationally more intensive.

5. Running the Code: Once config.yaml is set up, run the script. It will process the satellite observations and the model data according to the specified settings.

What the Code Does:

- Reads satellite observation data such as time, latitude, longitude, significant wave height, wind speed, and more from the given satellite .nc file.
- For each model data file in the specified directory, it will:
  - Create a regular grid based on the model's latitude and longitude.
  - Interpolate the model variables to the observation locations using the specified method.
  - Append the interpolated model data to the output file for each time step in the model data.
- The script maintains two time variables:
  - `time_sat`: Satellite observation time from the satellite data.
  - `time`: Time steps from the model data used for interpolation.
- Writes the combined data to the output netCDF file.

-author : Ghazal(Maryam) Mohammadpour 
email: ghazal.mohammadpour@noaa.gov

"""


import yaml
import netCDF4 as nc
import numpy as np
import os
from scipy.interpolate import LinearNDInterpolator
from pyresample import geometry, kd_tree

# Helper function to create the grid
def create_model_grid(model_lat, model_lon):
    lon, lat = np.meshgrid(model_lon, model_lat)
    return lat, lon

# Interpolate model data to observation locations
def interpolate_model_to_observation(obs_lat, obs_lon, model_grid, model_data, method):
    model_lat_grid, model_lon_grid = model_grid
    interpolated_data = {}

    if method == 'linear':
        model_points = np.column_stack((model_lat_grid.ravel(), model_lon_grid.ravel()))

        for var_data, var_name in model_data:
            var_values = var_data.ravel()
            interp = LinearNDInterpolator(model_points, var_values, fill_value=np.nan)
            interpolated_data[var_name] = interp(obs_lat, obs_lon)

    elif method == 'nearest':
        swath_def = geometry.SwathDefinition(lons=obs_lon, lats=obs_lat)
        grid_def = geometry.GridDefinition(lons=model_lon_grid, lats=model_lat_grid)

        for var_data, var_name in model_data:
            assert var_data.shape == model_lat_grid.shape, f"Data shape mismatch: var_data {var_data.shape}, model_lat_grid {model_lat_grid.shape}"
            resampled_data = kd_tree.resample_nearest(grid_def, var_data, swath_def, radius_of_influence=50000, fill_value=np.nan)
            interpolated_data[var_name] = resampled_data

    else:
        raise ValueError(f"Unknown interpolation method: {method}")

    return interpolated_data

if __name__ == "__main__":
    # Load configuration
    with open("config.yaml", 'r') as yamlfile:
        cfg = yaml.safe_load(yamlfile)

    obs_file = cfg["input"]["obs_file"]
    model_directory = cfg["input"]["model_directory"]
    output_file = cfg["output"]["file"]
    interp_method = cfg["interpolation"]["method"]

    # Read data from the input .nc file (satellite observations)
    with nc.Dataset(obs_file, 'r') as src:
        time_obs = src.variables['time'][:]
        latitude_obs = src.variables['latitude'][:]
        longitude_obs = src.variables['longitude'][:]
        sat_time = src.variables['sat_time'][:]
        hs = src.variables['hs'][:]
        hs_cal = src.variables['hs_cal'][:]
        wsp = src.variables['wsp'][:]
        wsp_cal = src.variables['wsp_cal'][:]
        sat_name = src.variables['sat_name'][:]

    # Write data to the output .nc file
    with nc.Dataset(output_file, 'w', format='NETCDF4') as dst:
        # Create dimensions
        dst.createDimension('time_obs', len(time_obs))
        dst.createDimension('obs', len(latitude_obs))
        dst.createDimension('sname', 1)

        # Create variables for satellite data
        time_var_sat = dst.createVariable('time_sat', 'f8', ('time_obs',))
        lat_var_sat = dst.createVariable('latitude_sat', 'f4', ('obs',))
        lon_var_sat = dst.createVariable('longitude_sat', 'f4', ('obs',))
        sat_time_var = dst.createVariable('sat_time_sat', 'f8', ('obs',))
        hs_var = dst.createVariable('hs_sat', 'f4', ('obs',))
        hs_cal_var = dst.createVariable('hs_cal_sat', 'f4', ('obs',))
        wsp_var = dst.createVariable('wsp_sat', 'f4', ('obs',))
        wsp_cal_var = dst.createVariable('wsp_cal_sat', 'f4', ('obs',))
        sat_name_var = dst.createVariable('sat_name_sat', str, ('sname',))

        # Assign attributes to satellite variables
        lat_var_sat.units = "degrees_north"
        lon_var_sat.units = "degrees_east"
        time_var_sat.units = "seconds since 1970-01-01 00:00:00"
        sat_time_var.units = "seconds since 1970-01-01 00:00:00"
        hs_var.units = "m"
        hs_cal_var.units = "m"
        wsp_var.units = "m/s"
        wsp_cal_var.units = "m/s"

        # Write satellite data
        lat_var_sat[:] = latitude_obs
        lon_var_sat[:] = longitude_obs
        time_var_sat[:] = time_obs
        sat_time_var[:] = sat_time
        hs_var[:] = hs
        hs_cal_var[:] = hs_cal
        wsp_var[:] = wsp
        wsp_cal_var[:] = wsp_cal
        sat_name_var[:] = sat_name

        # Create dimension and variable for model time
        dst.createDimension('time', None)
        time_var = dst.createVariable('time', 'f8', ('time',))

        # Initialize a list to collect model time data
        times = []

        # Process each model file for interpolation
        for model_file in sorted(os.listdir(model_directory)):
            if model_file.endswith(".nc"):
                with nc.Dataset(os.path.join(model_directory, model_file), 'r') as model_nc:
                    model_lat = model_nc.variables['latitude'][:]
                    model_lon = model_nc.variables['longitude'][:]
                    model_times = model_nc.variables['time'][:]

                    for t_index, model_time in enumerate(model_times):
                        times.append(model_time)  # Store model time
                        model_data = [
                            (model_nc.variables['HTSGW_surface'][t_index, :, :], 'HTSGW_surface'),
                            (model_nc.variables['WIND_surface'][t_index, :, :], 'WIND_surface')
                            # Add more variables if needed
                        ]

                        # Create model grid for interpolation
                        model_mesh_lat, model_mesh_lon = create_model_grid(model_lat, model_lon)

                        # Interpolate model data to observation locations
                        interpolated_model_data = interpolate_model_to_observation(
                            latitude_obs, longitude_obs, (model_mesh_lat, model_mesh_lon), model_data, interp_method)

                        # Write data for each variable and time step
                        for model_var_name, interpolated_values in interpolated_model_data.items():
                            # Define variable if it's not already defined
                            if model_var_name not in dst.variables:
                                model_var = dst.createVariable(
                                    model_var_name, 'f4', ('time', 'obs',), zlib=True)
                            else:
                                model_var = dst.variables[model_var_name]

                            # Write data for the current time step
                            model_var[t_index, :] = interpolated_values

        # Assign model time variable after all time steps have been processed
        time_var[:] = np.array(times)

    print(f"Data has been written to {output_file}")


