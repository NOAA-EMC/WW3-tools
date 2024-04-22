import netCDF4 as nc
import numpy as np
import os
import json
import netCDF4 as nc
import numpy as np
import glob
import os
import datetime

# Load configuration from JSON file
with open('/work/noaa/marine/Ghazal.Mohammadpour/Ghazal.Mohammadpour1/hr3_eva/WW3-tools/ww3tools/config.json') as f:
    config = json.load(f)

# Iterate over each model configuration
for model_config in config['models']:
    model_dir = model_config['model_dir']
    output_file = model_config['output_file']

    # Get a list of all .nc files in the model directory
    nc_files_model = glob.glob(os.path.join(model_dir, '*.nc'))

    # Initialize lists to store data
    all_model_hs = []
    all_obs_hs = []
    all_time = []
    all_model_wnd = []
    all_obs_wnd = []

    # Iterate over each file in the model directory
    for idx, file_path_model in enumerate(nc_files_model):
        print(f"Processing file {idx + 1} of {len(nc_files_model)}: {file_path_model}")
        # Load the NetCDF file
        dataset_model = nc.Dataset(file_path_model, 'r')

        try:
            # Extract necessary variables
            model_hs = dataset_model.variables['model_hs'][:]
            obs_hs = dataset_model.variables['obs_hs'][:]
            time = dataset_model.variables['time'][:]
            initial_condition_time = dataset_model.getncattr('initial_condition_time')
            model_wnd = dataset_model.variables['model_wnd'][:]
            obs_wnd = dataset_model.variables['obs_wnd'][:]

            # Close the NetCDF dataset
            dataset_model.close()

            # Calculate the desired time range (first 24 hours)
            desired_start_time = initial_condition_time + (120 * 3600)  # 48 hours in seconds
            desired_end_time = initial_condition_time + (144 * 3600)  # 72 hours in seconds

            # Find the indices of the time array corresponding to the desired time range
            indices = np.where((time >= desired_start_time) & (time <= desired_end_time))[0]

            # Extract the data for the desired time range
            model_hs_filtered = model_hs[indices]
            obs_hs_filtered = obs_hs[indices]
            time_filtered = time[indices]
            model_wnd_filtered = model_wnd[indices]
            obs_wnd_filtered = obs_wnd[indices]

            # Append data to the lists
            all_model_hs.append(model_hs_filtered)
            all_obs_hs.append(obs_hs_filtered)
            all_time.append(time_filtered)
            all_model_wnd.append(model_wnd_filtered)
            all_obs_wnd.append(obs_wnd_filtered)
        except AttributeError:
            print(f"Warning: 'initial_condition_time' attribute not found in file: {file_path_model}")

    # Concatenate data from all files
    all_model_hs = np.concatenate(all_model_hs)
    all_obs_hs = np.concatenate(all_obs_hs)
    all_time = np.concatenate(all_time)
    all_model_wnd = np.concatenate(all_model_wnd)
    all_obs_wnd = np.concatenate(all_obs_wnd)
    
    # Combine time, model_hs, and obs_hs into a single array for sorting
    combined_data = np.array(list(zip(all_time, all_model_hs, all_obs_hs, all_model_wnd, all_obs_wnd)))

    # Sort the combined data array based on time
    sorted_data = sorted(combined_data, key=lambda x: x[0])

    # Unpack the sorted data into separate arrays
    sorted_time, sorted_model_hs, sorted_obs_hs, sorted_model_wnd, sorted_obs_wnd = zip(*sorted_data)

    # Save filtered data to a new NetCDF file
    with nc.Dataset(output_file, 'w') as nc_out:
        nc_out.createDimension('time', len(all_time))
        nc_out.createVariable('time', 'f8', ('time',))
        nc_out.createVariable('model_hs', 'f8', ('time',))
        nc_out.createVariable('obs_hs', 'f8', ('time',))
        nc_out.createVariable('model_wnd', 'f8', ('time',))
        nc_out.createVariable('obs_wnd', 'f8', ('time',))
        nc_out.variables['time'][:] = sorted_time
        nc_out.variables['model_hs'][:] = sorted_model_hs
        nc_out.variables['obs_hs'][:] = sorted_obs_hs
        nc_out.variables['model_wnd'][:] = sorted_model_wnd
        nc_out.variables['obs_wnd'][:] = sorted_obs_wnd

