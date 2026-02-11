    
# The DEM file used below can be found at:	
# https://github.com/dengwirda/dem/releases/tag/v0.1.1
import os
import argparse
import configparser
import numpy as np
import netCDF4 as nc

import jigsawpy

from scipy.interpolate import RegularGridInterpolator
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

#from spacing import *

def parse_input_args():
    parser = argparse.ArgumentParser(description='Create a mask file with multiple methods.')
    parser.add_argument('--config', type=str, required=True, help='Path to the configuration file.')
    args = parser.parse_args()
    return args

def load_configuration(config_path):
    config = configparser.ConfigParser()
    config.read(config_path)

        # Creating a dictionary and populating it with configuration settings
    configurations = {
        'mesh_file': config.get('MeshSettings', 'mesh_file', fallback=''),
        'ww3_mesh_file':config.get('MeshSettings', 'WW3_mesh_file', fallback=''),
        'hfun_hmax': float(config.get('MeshSettings', 'hfun_hmax', fallback='100')),
        'black_sea': config.getint('CommandLineArgs', 'black_sea', fallback=3),
        'mask_file': config.get('CommandLineArgs', 'mask_file', fallback=''),
        'hmax': float(config.get('Spacing', 'hmax', fallback='100.0')),
        'hshr': float(config.get('Spacing', 'hshr', fallback='100')),
        'nwav': int(config.get('Spacing', 'nwav', fallback='400')),
        'hmin': float(config.get('Spacing', 'hmin', fallback='100.0')),
        'dhdx': float(config.get('Spacing', 'dhdx', fallback='0.05')),
        'dem_file': config.get('DataFiles', 'dem_file', fallback='')
    }
    return configurations

ISOLATED = 30000.  # min surface area [km^2]

# just global objects, to keep things simple...
geom = jigsawpy.jigsaw_msh_t()
spac = jigsawpy.jigsaw_msh_t()
mesh = jigsawpy.jigsaw_msh_t()
opts = jigsawpy.jigsaw_jig_t()



def inject_dem():

#    args = parse_input_args()
#    configurations = load_configuration(args.config)

#-- remap a DEM on to the vertices of the mesh

    print("*inject-dem...")

        # Load the DEM file from the config
    #dem_file = configurations['dem_file']
    #data = nc.Dataset(dem_file,"r")

#    data = nc.Dataset("../RWPS/Data/RTopo_2_0_4_GEBCO_v2023_60sec_pixel.nc","r")
#    data = nc.Dataset("/scratch3/NCEPDEV/climate/Keston.Smith/RWPS/Data/RTopo_2_0_4_GEBCO_v2023_60sec_pixel.nc","r")
    data = nc.Dataset("../Data/Bathymetry/RTopo_2_0_4_GEBCO_v2024_60sec_pixel.nc","r")
    xlon = np.asarray(data["lon"][:])
    ylat = np.asarray(data["lat"][:])
    elev = np.asarray(data["bed_elevation"][:]) + \
           np.asarray(data["ice_thickness"][:])
        
    xmid = 0.5 * (xlon[:-1:] + xlon[1::])
    ymid = 0.5 * (ylat[:-1:] + ylat[1::])
        
    ffun = RegularGridInterpolator(
        (ymid, xmid), elev, 
        bounds_error=False, fill_value=None)

    vert = mesh.point["coord"]
    
    mids =(vert[mesh.tria3["index"][:, 0], :] +
           vert[mesh.tria3["index"][:, 1], :] +
           vert[mesh.tria3["index"][:, 2], :] 
          ) / 3.0 

    vsph = jigsawpy.R3toS2(geom.radii, vert)
    vsph*= 180. / np.pi
    
    mesh.value = ffun((vsph[:, 1], vsph[:, 0]))
    
    msph = jigsawpy.R3toS2(geom.radii, mids)
    msph*= 180. / np.pi

    mesh.vmids = ffun((msph[:, 1], msph[:, 0]))
    
    # save lon-lat at cell centres
    mesh.smids = np.zeros(
        (mesh.tria3.size, 2), dtype=np.float64)
    mesh.smids[:, 0] = msph[:, 0]
    mesh.smids[:, 1] = msph[:, 1]


def tri_to_tri(tria):

#-- return tria-to-tria adj. as a sparse graph

    # non-unique edges in tris
    edge = np.empty((0, 2), dtype=np.int32)
    tris = np.empty((0), dtype=np.int32)
    edge = np.concatenate((edge, 
        tria[:, (0, 1)]), axis=0)
    tris = np.concatenate((tris, 
        np.arange(0, tria.shape[0])))
        
    edge = np.concatenate((edge, 
        tria[:, (1, 2)]), axis=0)
    tris = np.concatenate((tris, 
        np.arange(0, tria.shape[0])))
        
    edge = np.concatenate((edge, 
        tria[:, (2, 0)]), axis=0)
    tris = np.concatenate((tris, 
        np.arange(0, tria.shape[0])))
        
    # which edges match to which?
    edge = np.sort(edge, axis=1)
    imap = np.argsort(edge[:, 1], kind="stable")
    edge = edge[imap, :]
    tris = tris[imap]
    imap = np.argsort(edge[:, 0], kind="stable")
    edge = edge[imap, :]
    tris = tris[imap]
    
    diff = edge[1::, :] - edge[:-1:, :]

    same = np.argwhere(np.logical_and.reduce((
        diff[:, 0] == 0, 
        diff[:, 1] == 0))).ravel()
        
    # tris[same] and tris[same+1] share
    rows = np.concatenate((
        tris[same], tris[same+1]))
    cols = np.concatenate((
        tris[same+1], tris[same]))
    data = np.ones(rows.size, dtype=np.int8)
    
    # ith tri is adj. to tri in ith row
    return csr_matrix((data, (rows, cols)))


def filter_dry(mesh, mask):

#-- require dry cells > 1 dry edge, and large

    print("*filter-dry...")

    filt = np.logical_not(mask)
    tris = np.argwhere(filt).ravel()

    conn = tri_to_tri(mesh.tria3["index"][filt, :])
        
    # require dry to be adj. >=1 dry cell
    isol = np.sum(conn, axis=1) <= 1
    isol = np.ravel(isol)
#KWS COMMENT BELOW
    
    # delete groups of dry if too small
    nprt, part = connected_components(
        conn, directed=False, return_labels=True)

    tris = np.argwhere(filt).ravel()
    for iprt in range(nprt):
        itri = np.argwhere(part == iprt)
        if (itri.size <= 2): mask[tris[itri]] = True
    
#KWS COMMENT ABOVE
 
    # otherwise mark isolated cell as ocn
    mask[tris[isol]] = True

    return mask


def filter_wet(mesh, mask):

#-- require wet cells > 1 wet edge, and large

    print("*filter-wet...")

    tris = np.argwhere(mask).ravel()

    conn = tri_to_tri(mesh.tria3["index"][mask, :])

    # require wet to be adj. >=1 wet cell
    isol = np.sum(conn, axis=1) <= 1
    isol = np.ravel(isol)
    
    # delete groups of wet if too small
    nprt, part = connected_components(
        conn, directed=False, return_labels=True)

    area = jigsawpy.trivol2(
        mesh.point["coord"], 
        mesh.tria3["index"][mask, :])

    for iprt in range(nprt):
        itri = np.argwhere(part == iprt)
        asum = np.sum(area[itri])
       #print(asum)
        if asum < ISOLATED: mask[tris[itri]] = False
    
    # otherwise mark isolated cell as dry
    mask[tris[isol]] = False
    
    return mask


def filter_ocn():

#    args = parse_input_args()
#    configurations = load_configuration(args.config)

#-- use the remapped elev. to keep ocean cells

    print("*filter-ocn...")
#KWS not sure what this is
#    elev =(mesh.value[mesh.tria3["index"][:, 0]]
#         + mesh.value[mesh.tria3["index"][:, 1]]
#         + mesh.value[mesh.tria3["index"][:, 2]]
#         + mesh.vmids) / 4.0
    elev =(mesh.value[mesh.tria3["index"][:, 0]]
         + mesh.value[mesh.tria3["index"][:, 1]]
         + mesh.value[mesh.tria3["index"][:, 2]]
         ) / 3.0
    # Define the Caspian Sea region
    caspian_lat_min = 34.5
    caspian_lat_max = 50.0
    caspian_lon_min = 44.5
    caspian_lon_max = 55.5
    
    # Define the black Sea region
    blacksea_lat_min = 40 
    blacksea_lat_max = 47.25
    blacksea_lon_min = 26.15
    blacksea_lon_max = 41.5
    
   # Define the additional region1
    additional_lat_min = 39.95
    additional_lat_max = 40.6
    additional_lon_min = 26
    additional_lon_max = 26.8

   # Define the additional region2
    additional_lat_min2 = 40.3
    additional_lat_max2 = 40.6
    additional_lon_min2 = 26.8
    additional_lon_max2 = 30

   # Define the additional region3
    additional_lat_min3 = 40.6
    additional_lat_max3 = 41.25
    additional_lon_min3 = 28.9
    additional_lon_max3 = 29.1
    
    """
# Print the elevations in the additional region
    additional_elev = elev[np.logical_and.reduce((
        mesh.smids[:, 1] >= additional_lat_min,
        mesh.smids[:, 1] <= additional_lat_max,
        mesh.smids[:, 0] >= additional_lon_min,
        mesh.smids[:, 0] <= additional_lon_max
    ))]
    print("Elevations in the additional region:")
    print(additional_elev)

# Print the elevations in the additional region2
    additional_elev2 = elev[np.logical_and.reduce((
        mesh.smids[:, 1] >= additional_lat_min2,
        mesh.smids[:, 1] <= additional_lat_max2,
        mesh.smids[:, 0] >= additional_lon_min2,
        mesh.smids[:, 0] <= additional_lon_max2
    ))]
    print("Elevations in the additional region2:")
    print(additional_elev2)

# Print the elevations in the additional region3
    additional_elev3 = elev[np.logical_and.reduce((
        mesh.smids[:, 1] >= additional_lat_min3,
        mesh.smids[:, 1] <= additional_lat_max3,
        mesh.smids[:, 0] >= additional_lon_min3,
        mesh.smids[:, 0] <= additional_lon_max3
    ))]
    print("Elevations in the additional region3:")
    print(additional_elev3)

    """

    # zssh, to cull elev. against
    surf = np.zeros(elev.shape, dtype=np.float32)
    # Update the surf array to include both regions
    # Define the Caspian Sea region
    caspian_region = np.logical_and.reduce((
        mesh.smids[:, 1] >= caspian_lat_min,
        mesh.smids[:, 1] <= caspian_lat_max,
        mesh.smids[:, 0] >= caspian_lon_min,
        mesh.smids[:, 0] <= caspian_lon_max
    ))
    
    blacksea_region = np.logical_and.reduce((
        mesh.smids[:, 1] >= blacksea_lat_min,
        mesh.smids[:, 1] <= blacksea_lat_max,
        mesh.smids[:, 0] >= blacksea_lon_min,
        mesh.smids[:, 0] <= blacksea_lon_max
    ))

    keep = elev <= surf  # only keep tri with wet elev
   
   
   
    # iterate on dry cells until none "isolated" 
    #KWS remove filter dry for Non Global mesh
#   if (False):
#   	knum = np.count_nonzero(keep)
#	while (True):   
#		keep = filter_dry(mesh, keep)
#        	if (np.count_nonzero(keep) == knum): break
#		knum = np.count_nonzero(keep)
#	mesh.tria3 = mesh.tria3[keep]

    # iterate on wet cells until none "isolated"
    keep = np.ones(mesh.tria3.size, dtype=bool)
    knum = np.count_nonzero(keep)
    while (True):
        keep = filter_wet(mesh, keep)
        if (np.count_nonzero(keep) == knum): break
        knum = np.count_nonzero(keep)
    
    mesh.tria3 = mesh.tria3[keep]

    # delete unused vertices and reindex
    ifwd = np.unique(mesh.tria3["index"].ravel())
    
    irev = np.zeros(mesh.point.size, dtype=np.int32)
    irev[ifwd] = np.arange(ifwd.size, dtype=np.int32)

    mesh.point = mesh.point[ifwd]
    mesh.value = mesh.value[ifwd]
    mesh.tria3["index"] = irev[mesh.tria3["index"]]
    return mesh

def write_gmsh_mesh(filename, node_data, tri):
    num_nodes = node_data.shape[0]

    # Write the mesh to the file
    with open(filename, 'w') as fileID:
        fileID.write("$MeshFormat\n")
        fileID.write("2 0 8\n")
        fileID.write("$EndMeshFormat\n")
        fileID.write("$Nodes\n")
        fileID.write(str(num_nodes) + "\n")

        for i in range(num_nodes):
            fileID.write(
                f"{i + 1}  {node_data[i, 0]:5.5f} {node_data[i, 1]:5.5f} {node_data[i, 2]:5.5f}\n"
            )

        fileID.write("$EndNodes\n")
        fileID.write("$Elements\n")
        num_elements = len(tri)
        fileID.write(str(num_elements) + "\n")

        m = 0
        for i in range(len(tri)):
            m += 1
            fileID.write(f"{m} 2 3 0 {i+1} 0 {tri[i][0]} {tri[i][1]} {tri[i][2]}\n")

        fileID.write("$EndElements\n")

def interp_bathymetry_v0():


#-- remap a DEM on to the vertices of the mesh

    print("*interp-batyhmetry...")

        # Load the DEM file from the config
    dem_file = "../RWPS/Data/RTopo_2_0_4_GEBCO_v2023_60sec_pixel.nc"

    data = nc.Dataset(dem_file,"r")

    xlon = np.asarray(data["lon"][:])
    ylat = np.asarray(data["lat"][:])
    elev = np.asarray(data["bed_elevation"][:]) + \
           np.asarray(data["ice_thickness"][:])
        
    xmid = 0.5 * (xlon[:-1:] + xlon[1::])
    ymid = 0.5 * (ylat[:-1:] + ylat[1::])
        
    ffun = RegularGridInterpolator((ymid, xmid), elev,bounds_error=False, fill_value=None)

    vert = mesh.point["coord"]
    hdown = ffun( (vert[:, 1], vert[:, 0]) )

    vert[:,2]=hdown
    mesh.point["coord"][:,2]=vert[:,2] 
    return mesh


def interp_bat():

#-- remap a DEM on to the vertices of the mesh

    print("*interp-batyhmetry...")
        
    ffun = RegularGridInterpolator((topo.ygrid, topo.xgrid),topo.value,bounds_error=False, fill_value=None)
    vert = mesh.point["coord"]
    hdown = ffun( (vert[:, 1], vert[:, 0]) )
    mesh.value=hdown
    return mesh

def do_nothing( ):
    print("the end")

