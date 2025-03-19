""" Jigsaw meshes for WW3 with global bathymetry
"""

# Authors: Ali Salimi-Tarazouj, Darren Engwirda

# This script will create mesh spacing based on user defined windows in json format or based on polygons in shapefile format

import configparser
import numpy as np
import netCDF4 as nc
import argparse
import json
import geopandas as gpd
from shapely.geometry import Point

def parse_input_args():
    parser = argparse.ArgumentParser(description='Create a mask file with multiple methods.')
    parser.add_argument('--config', type=str, required=True, help='Path to the configuration file.')
    args = parser.parse_args()
    return args

def load_configuration(config_path):
    config = configparser.ConfigParser()
    config.read(config_path)
    try:
        windows = json.loads(config['DataFiles']['window_file']) if 'window_file' in config['DataFiles'] else None
    except json.JSONDecodeError:
        print("Error parsing windows data")
        windows = None

    try:
        shapefiles = json.loads(config['DataFiles']['shape_file']) if 'shape_file' in config['DataFiles'] else None
    except json.JSONDecodeError:
        print("Error parsing shapefiles data")
        shapefiles = None

    if shapefiles:
        shapefiles = [(shp['path'], shp['scale']) for shp in shapefiles]

    dem_file = config['DataFiles']['dem_file'] if 'dem_file' in config['DataFiles'] else None

    scaling_settings = {key: float(value) for key, value in config['ScalingSettings'].items()} if 'ScalingSettings' in config else {}

    return windows, shapefiles, dem_file, scaling_settings

def create_mask_file(data_filename, output_filename, windows=None, shapefiles=None, default_scale=20):
    # Load DEM data
    data = nc.Dataset(data_filename, "r")
    xlon = np.asarray(data["lon"][:])
    ylat = np.asarray(data["lat"][:])
    elev = np.asarray(data["bed_elevation"][:], dtype=np.float32) + np.asarray(data["ice_thickness"][:], dtype=np.float32)

    # Compute midpoints for longitude and latitude
    xmid = 0.5 * (xlon[:-1] + xlon[1:])
    ymid = 0.5 * (ylat[:-1] + ylat[1:])
    xmat, ymat = np.meshgrid(xmid, ymid)

    # Scaling settings from the config file
    upper_bound = scaling_settings.get('upper_bound')
    middle_bound = scaling_settings.get('middle_bound')
    lower_bound = scaling_settings.get('lower_bound')
    scale_north = scaling_settings.get('scale_north')
    scale_middle = scaling_settings.get('scale_middle')
    scale_south_upper = scaling_settings.get('scale_south_upper')
    scale_south_lower = scaling_settings.get('scale_south_lower')


    # Apply scaling based on latitude bands
    scal = np.where(ymat > upper_bound, scale_north,
                    np.where(ymat > middle_bound, scale_middle,
                             scale_south_upper + (scale_south_lower - scale_south_upper) *
                             (ymat - middle_bound) / (lower_bound - middle_bound)))

    # Apply window-based refinement if provided
    if windows:
        for window in windows:
            window_mask_lon = np.logical_and(xmat >= window["min_lon"], xmat <= window["max_lon"])
            window_mask_lat = np.logical_and(ymat >= window["min_lat"], ymat <= window["max_lat"])
            window_mask = np.logical_and(window_mask_lon, window_mask_lat)

            mask_elevation = np.logical_and(elev >= -4., elev <= +8.)
            mask_final = np.logical_and(mask_elevation, window_mask)
            scal[mask_final] = window["hshr"]

    # Process shapefiles if provided
    if shapefiles:
        for shapefile, scale in shapefiles:
            gdf = gpd.read_file(shapefile)
            minx, miny, maxx, maxy = gdf.total_bounds
            # Subset the meshgrid
            mask_x = (xmid >= minx) & (xmid <= maxx)
            mask_y = (ymid >= miny) & (ymid <= maxy)
            x_sub, y_sub = xmid[mask_x], ymid[mask_y]
            xmat_sub, ymat_sub = np.meshgrid(x_sub, y_sub)

            # Create GeoDataFrame for the points in the subset
            points_sub = [Point(x, y) for x, y in zip(xmat_sub.flatten(), ymat_sub.flatten())]
            points_gdf_sub = gpd.GeoDataFrame(geometry=points_sub)
            points_gdf_sub.set_crs(gdf.crs, inplace=True)

            # Spatial join only for the subset area
            overlay_sub = gpd.sjoin(points_gdf_sub, gdf, how="left", op='intersects')
            overlay_sub['scale'] = overlay_sub['index_right'].apply(lambda x: scale if x >= 0 else default_scale)
            scales_sub = overlay_sub['scale'].to_numpy().reshape(xmat_sub.shape)

            # Map the scales back to the original full matrix
            x_indices = np.searchsorted(xmid, x_sub)
            y_indices = np.searchsorted(ymid, y_sub)
            scal[np.ix_(mask_y, mask_x)] = scales_sub

    # Write results to a NetCDF file
    data_out = nc.Dataset(output_filename, "w")
    data_out.createDimension("nlon", xmid.size)
    data_out.createDimension("nlat", ymid.size)
    var = data_out.createVariable("val", "f4", ("nlat", "nlon"))
    var[:, :] = scal
    data_out.close()

if __name__ == "__main__":
    args = parse_input_args()
    windows, shapefiles, dem_file, scaling_settings = load_configuration(args.config)
    create_mask_file(dem_file, "wmask.nc", windows, shapefiles)

