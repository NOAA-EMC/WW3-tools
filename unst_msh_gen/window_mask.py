""" Jigsaw meshes for WW3 with global bathymetry
"""

# Authors: Ali Salimi, Darren

# This script will create mesh spacing based on user defined windows in json format.

import numpy as np
import netCDF4 as nc
import argparse
import json

def create_mask_file(windows):
    # Open the data file
    data = nc.Dataset("RTopo_2_0_4_GEBCO_v2023_60sec_pixel.nc", "r")

    # Extract longitude and latitude arrays
    xlon = np.asarray(data["lon"][:])
    ylat = np.asarray(data["lat"][:])
    elev = np.asarray(data["bed_elevation"][:], dtype=np.float32) + np.asarray(data["ice_thickness"][:], dtype=np.float32)

    # Compute midpoints for longitude and latitude
    xmid = 0.5 * (xlon[:-1] + xlon[+1:])
    ymid = 0.5 * (ylat[:-1] + ylat[+1:])

    # Create a meshgrid of midpoints
    xmat, ymat = np.meshgrid(xmid, ymid)

    # Initialize scal with a default condition
    #scal = np.where(ymat > 50, 9, np.where(ymat < -20, 30, 20))

    # Define the boundaries and scaling values: The globe is divided to three regions based on latitude and corresponding resolutions
    upper_bound = 50       # Lattitude that starts the upper section (ie lat > upper_bound) 
    middle_bound = -20     # Determines the boundary between the upper  and middle sections.  The middle section is for upper_bound > lat > middle_bound 
    lower_bound = -90      # Boundary of lowest section, likely -90,   Lower section is:   middle_bound > lat > lower_bound
    scale_north = 9        # mesh resolution  in km for upper section ( lat > upper_bound )
    scale_middle = 20      # mesh resolution in km for  middle section for middle_bound < lat < upper_boun
    # Mesh resolution of the lower section is linear from scale_south_upper to scale_south_upper
    scale_south_upper = 30 # km resolution for upper south/lower section 
    scale_south_lower = 9  # km mesh resolution at lower south/lower section 

    # Calculate the scaling using conditions
    scal = np.where(ymat > upper_bound,
                    scale_north,  # Use 9 km for ymat > 50
                    np.where(ymat > middle_bound,
                             scale_middle,  # Use 20 km for -20 < ymat <= 50
                             # Gradual decrease from 30 km to 9 km as latitude decreases from -20 to -90
                             scale_south_upper + (scale_south_lower - scale_south_upper) * 
                             (ymat - middle_bound) / (lower_bound - middle_bound)
                            ))

    # Process each window
    for window in windows:
        # Apply the 'hshr' value for the current window where conditions are met
        window_mask_lon = np.logical_and(xmat >= window["min_lon"], xmat <= window["max_lon"])
        window_mask_lat = np.logical_and(ymat >= window["min_lat"], ymat <= window["max_lat"])
        window_mask = np.logical_and(window_mask_lon, window_mask_lat)
        
        mask_elevation = np.logical_and(elev >= -4., elev <= +8.)
        mask_final = np.logical_and(mask_elevation, window_mask)
        scal[mask_final] = window["hshr"]

    # Create a new NetCDF file to store the mask
    data_out = nc.Dataset("wmask.nc", "w")
    data_out.createDimension("nlon", xmid.size)
    data_out.createDimension("nlat", ymid.size)
    if "val" not in data_out.variables.keys():
        data_out.createVariable("val", "f4", ("nlat", "nlon"))
    data_out["val"][:, :] = scal
    data_out.close()

def parse_windows_from_args():
    parser = argparse.ArgumentParser(description='Create a mask file with multiple windows.')
    parser.add_argument('--windows', type=str, required=True,
                        help='JSON string with window definitions. Example: \'[{"min_lon": -130, "max_lon": -64, "min_lat": 24.5, "max_lat": 47.5, "hshr": 5}]\'')
    
    args = parser.parse_args()
    windows = json.loads(args.windows)
    return windows

if __name__ == "__main__":
    windows = parse_windows_from_args()
    create_mask_file(windows)

