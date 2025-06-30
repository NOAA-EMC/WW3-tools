#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 15:17:57 2024

@author: Ali Salimi-Tarazouj

This script originally was developed by Steven Brus (@sbrus89) and revised by Ali Salimi-Tarazouj to work efficiently for very high resolution meshes
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.ticker as mticker
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import argparse


def read_gmsh(filename):
    #purpose: this function reads a gmsh file and returns node and element information
    # additional information about the gmsh file format can be found
    # in many places including here: http://gmsh.info/dev/doc/texinfo/gmsh.pdf
    # The gmsh format is what the WW3 model uses to define unstructured grids

    #input:
    # filename - name of gmsh file

    #output:
    #xy    -  x/y or lon/lat of nodes
    #depth - depth value at node points
    #ect   - element connection table
    #bnd   - list of boundary nodes

    with open(filename, 'r') as f:
        # Skip mesh format lines
        for _ in range(4):  # Skip 4 lines directly (including $Nodes)
            next(f)
        
        # Read number of nodes
        nn = int(next(f).strip())
        
        # Initialize node coordinate and depth arrays
        xy = np.zeros((nn, 2), dtype=np.double)
        depth = np.zeros(nn, dtype=np.double)
        
        # Read coordinates and depths
        for i in range(nn):
            line = next(f).split()
            idx = int(line[0]) - 1
            x_coord = float(line[1])
            xy[idx, 0] = x_coord - 360 if x_coord > 180 else x_coord
            xy[idx, 1] = float(line[2])
            depth[idx] = float(line[3])
        
        # Skip '$EndNodes' and '$Elements'
        next(f)
        next(f)
        
        # Read number of elements
        ne = int(next(f).strip())
        
        # Initialize temporary arrays to read in element info
        ecttemp = np.zeros((ne, 3), dtype=np.int32)
        bndtemp = []
        
        elem_count = 0
        for _ in range(ne):
            line = next(f).split()
            eltype = int(line[1])
            if eltype == 15:
                bndtemp.append(int(line[5]) - 1)
            else:
                ecttemp[elem_count, :] = [int(line[6]) - 1, int(line[7]) - 1, int(line[8]) - 1]
                elem_count += 1
        
        # Trim the ect array to the actual number of elements
        ect = ecttemp[:elem_count, :]
        bnd = np.array(bndtemp, dtype=np.int32)
    
    return xy, depth, ect, bnd

def calc_elm_size(xy, ect):

    #purpose: Calculate element size, by calculating the distance between each node on the triangle, then
   # saving the minimum or maximum distance for each node

   #input:
   # xy  -  x/y or lon/lat of nodes
   # ect - element connection table

   #output:
   # distmin(number of nodes)  - the minimum distance between this node and any connected node point
   # distmax(number of nodes)  - the maximum distance between this node and any connected node point

    nn = len(xy)
    ne = len(ect)
    radiusofearth = 6378.137  # radius of earth at equator
    d2r = np.pi / 180  # degrees to radians conversion factor

    # Initialize arrays to store min and max distances
    distmin = np.full(nn, np.inf)
    distmax = np.zeros(nn)

    # Precompute cosines and sines for latitudes
    cos_lat = np.cos(xy[:, 1] * d2r)
    sin_lat = np.sin(xy[:, 1] * d2r)

    # Calculate distances for all elements
    for i, j, k in ect:
        # Extract coordinates
        lon_i, lon_j, lon_k = xy[i, 0], xy[j, 0], xy[k, 0]
        cos_lat_i, cos_lat_j, cos_lat_k = cos_lat[i], cos_lat[j], cos_lat[k]
        sin_lat_i, sin_lat_j, sin_lat_k = sin_lat[i], sin_lat[j], sin_lat[k]

        # Calculate distances using the Haversine formula
        dist_ij = 2 * radiusofearth * np.arcsin(np.sqrt(np.sin((xy[j, 1] - xy[i, 1]) * d2r / 2) ** 2 +
                                                       cos_lat_i * cos_lat_j * np.sin((lon_j - lon_i) * d2r / 2) ** 2))
        dist_jk = 2 * radiusofearth * np.arcsin(np.sqrt(np.sin((xy[k, 1] - xy[j, 1]) * d2r / 2) ** 2 +
                                                       cos_lat_j * cos_lat_k * np.sin((lon_k - lon_j) * d2r / 2) ** 2))
        dist_ki = 2 * radiusofearth * np.arcsin(np.sqrt(np.sin((xy[i, 1] - xy[k, 1]) * d2r / 2) ** 2 +
                                                       cos_lat_k * cos_lat_i * np.sin((lon_i - lon_k) * d2r / 2) ** 2))

        # Update min and max distances for each node
        distmin[i], distmax[i] = min(distmin[i], dist_ij, dist_ki), max(distmax[i], dist_ij, dist_ki)
        distmin[j], distmax[j] = min(distmin[j], dist_ij, dist_jk), max(distmax[j], dist_ij, dist_jk)
        distmin[k], distmax[k] = min(distmin[k], dist_jk, dist_ki), max(distmax[k], dist_jk, dist_ki)

    return distmin, distmax

def create_mask(xy, ect):
    # Extract x-coordinates of the nodes for each element
    x_coords = xy[ect, 0]  # ect indexes rows in xy, column 0 is x-coordinates

    # Check for cross-dateline elements
    # Calculate signs and check if all are the same for each element
    signs = np.sign(x_coords)
    uniform_sign = (np.ptp(signs, axis=1) == 0)  # Change in peak-to-peak (max-min) across rows; if 0, all signs are the same

    # Check for elements near the dateline, i.e., abs(x) < 10
    near_dateline = (np.abs(x_coords) < 10).any(axis=1)  # Check any coordinate < 10 in absolute value for each element

    # Combine conditions: elements are valid if they do not cross the dateline and are not near the dateline
    mask = ~(uniform_sign | near_dateline)

    return mask.astype(int)  # Convert boolean mask to integer (1s and 0s)



def setup_plot(ax, extent=[-180, 180, -90, 90]):
    """ Common plot setup function """
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE, zorder=101)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlines = False
    gl.ylines = False
    gl.xlocator = mticker.FixedLocator(np.linspace(extent[0], extent[1], 7))
    gl.ylocator = mticker.FixedLocator(np.linspace(extent[2], extent[3], 7))
    gl.xformatter = LongitudeFormatter()
    gl.yformatter = LatitudeFormatter()

def plot_eleminfo(plotdescriptor, xy, ect, distmin, distmax, depth, highlight_nodes=None, named_points=None):

    print('Grid')
    print(plotdescriptor)
    print('Min/max of distmin:', np.min(distmin), np.max(distmin))
    print('Min/max of distmax:', np.min(distmax), np.max(distmax))
    print('Min/max of bathy:', np.min(depth), np.max(depth))

    mask = create_mask(xy,ect)
  
  #create one tringulation and use it for multiple plots: 
    triang=tri.Triangulation(xy[:,0],xy[:,1],triangles=ect, mask=mask)
    vpltmin=1
    vpltmax= 30
        

    # Shared figure setup
    figsize = [18.0, 9.0]
    plots = [
        ('elm', 'k-', 'Element outlines', None, 'jet'),
        ('maxsize', distmax, 'Element size in km', 'gouraud', 'jet'),
        ('minsize', distmin, 'Element size in km', 'gouraud', 'jet'),
        ('bathy', depth, 'Bathymetry in m', None, 'jet')
    ]

    for suffix, data, label, shading, cmap in plots:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        setup_plot(ax)
        if suffix == 'elm':
            cf = ax.triplot(triang, 'k-', linewidth=0.5)
        else:
            cf = ax.tripcolor(triang, data, cmap=cmap, shading=shading if shading else 'flat', vmin=vpltmin, vmax=vpltmax if data is not depth else None)
            plt.colorbar(mappable=cf, label=label)

        if highlight_nodes is not None and suffix == 'elm':
            ax.scatter(xy[highlight_nodes, 0], xy[highlight_nodes, 1], color='red', s=100, zorder=102, transform=ccrs.PlateCarree())
        # --- ADD THIS BLOCK FOR NAMED POINTS ---
        if named_points is not None and len(named_points) > 0:
            for lon, lat, name in named_points:
                if lon > 180 : 
                    lon = lon - 360 
                ax.scatter(lon, lat, color='black', s=10, marker='o', zorder=103, transform=ccrs.PlateCarree())
                ax.text(lon, lat, name, fontsize=10, color='black', zorder=104, transform=ccrs.PlateCarree(), ha='left', va='bottom')
        # ---------------------------------------


        plt.savefig(f'{suffix}_{plotdescriptor}_{name}.png')
        #plt.close()

def main():
    parser = argparse.ArgumentParser(description='Plot mesh information from a gmsh file.')
    parser.add_argument('--filename', type=str, required=True, help='Path to the gmsh file')
    parser.add_argument('--descriptor', type=str, required=True, help='Descriptor for output plot filenames')
    args = parser.parse_args()

    filename = args.filename
    descriptor = args.descriptor

    xy, depth, ect, bnd =  read_gmsh(filename)
    distmin, distmax = calc_elm_size(xy, ect)
    highlighted_nodes = []  # Replace with your actual node indices, e.g. [12776, 13923]
    named_points = [] # Replace with similar list as seen below: 
    #named_poins = [ 
    #(-165.475,   64.473,  "46265"),
    #(-157.750,   21.480,  "51207"),
    # # ...
    #]
    plot_eleminfo(descriptor, xy, ect, distmin, distmax, depth, highlight_nodes=highlighted_nodes, named_points=named_points)

if __name__ == "__main__":
    main()
