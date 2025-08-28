
""" Jigsaw meshes for WW3 with global bathymetry
"""

# Authors: Ali Salimi-Tarazouj, Darren Engwirda

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

from spacing import *

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

def create_msh():

#-- create a simple uniform mesh for the globe
    
    args = parse_input_args()
    configurations = load_configuration(args.config)

    print("*create-msh...")

    opts.geom_file = "geom.msh"  #saves the geometry info for jigsaw
    opts.hfun_file = "spac.msh"  #saves the final mesh spacing info 
    opts.jcfg_file = "opts.jig"  #jigsaw ctlr file
    
    geom.mshID = "ellipsoid-mesh"
    geom.radii = np.full(
        3, 6.371E+003, dtype=geom.REALS_t)

    jigsawpy.savemsh(opts.geom_file, geom)
    
    create_siz()
    
    jigsawpy.savemsh(opts.hfun_file, spac)
    
    # solve |dh/dx| constraints in spacing
    jigsawpy.cmd.marche(opts, spac)
    
    opts.mesh_file = configurations['mesh_file']  #jigsaw format mesh file
    
    opts.hfun_scal = "absolute"
    opts.hfun_hmax = configurations['hfun_hmax']           # global maximum mesh resolution (similar to hmax)
    opts.mesh_dims = +2             # 2-dim. simplexes
    opts.optm_iter = +64            # number of itereation for the optimization
    opts.optm_cost = "skew-cos"

    jigsawpy.cmd.jigsaw(opts, mesh)
    
    
def create_siz():

    args = parse_input_args()
    configurations = load_configuration(args.config)

    #-- create mesh spacing function for the globe: for uniform mesh hmax = hshr = hmin

    hmax = configurations['hmax'] # maximum spacing [km] 
    hshr = configurations['hshr']   # shoreline spacing
    nwav = configurations['nwav']   # number of cells per sqrt(g*H)
    hmin = configurations['hmin']  # minimum spacing
    dhdx = configurations['dhdx']  # allowable spacing gradient: for more gradual transition use lower value
    mask_file = configurations['mask_file'] #user defined scaling file
    # Load the DEM file from the config
    dem_file = configurations['dem_file']

    data = nc.Dataset(dem_file,"r")

    xlon = np.asarray(data["lon"][:])
    ylat = np.asarray(data["lat"][:])
    elev = np.asarray(data["bed_elevation"][:]) + \
           np.asarray(data["ice_thickness"][:])
           
    land = form_land_mask_connect(elev, edry=2) >= 1
    high = form_land_mask_connect(elev, edry=8) >= 1

#-- init. h(x) data: impose global "reachable" land mask

    hmat = np.full(
        (elev.shape[:]), hmax, dtype=spac.FLT32_t)
        
    hmat[land] = hmax

    if (nwav > 0.0):
        hmat = np.minimum(
            hmat, swe_wavelength_spacing(
                elev, land, nwav, hmin, hmax))
    
#-- final h(x) data: impose global "shoreline" min. val.
   
    hmat[high] = hmax
   
    hmat = setup_shoreline_pixels(hmat, land, hshr)
   
#-- apply user-defined scaling: multiply h(x) by mask array

    if mask_file:
        hmat = scale_spacing_via_mask(args, hmat)
        print("Scaling applied using mask_file:", mask_file)
    else:
    # Handle case where mask_file is not provided
        print("No mask file provided. Proceeding without scaling...")

#-- and a little nonlinear smoothing
    
    filt = filter_pixels_harmonic(hmat, exp=2)
    hmat = np.minimum(hmat, filt)
        
    filt = filter_pixels_harmonic(hmat, exp=1)
    hmat = np.minimum(hmat, filt)

    hmat = np.asarray(remap_pixels_to_corner(hmat), 
                      dtype=spac.FLT32_t)
    
#-- pack h(x) data to jigsaw datatype: average pixel-to-
#-- node, careful with periodic BCs.
    
    spac.mshID = "ellipsoid-grid"
    spac.radii = geom.radii
    spac.xgrid = xlon * np.pi / 180.
    spac.ygrid = ylat * np.pi / 180.

    xmat, ymat = np.meshgrid(
        spac.xgrid, spac.ygrid, sparse=True)

#-- keep high-res. only in a guassian-ish "zoom" region

    ymid = 41.5 * np.pi / 180.
    xmid = 30.5 * np.pi / 180.

    zoom = +100.0 - 99.0 * np.exp(-(
        6.75 * (xmat - xmid) ** 2 +
        12.5 * (ymat - ymid) ** 2) ** 2)
    
    spac.value = hmat*zoom 
    spac.slope = np.array(dhdx)
    spac.value = np.minimum(hmax, spac.value)
    
#-- save spacing to a netcdf, for viz. in e.g. paraview
    
    data = nc.Dataset("spac.nc", "w")
    data.createDimension("nlon", spac.xgrid.size)
    data.createDimension("nlat", spac.ygrid.size) 

    if ("val" not in data.variables.keys()):
        data.createVariable("val", "f4", ("nlat", "nlon"))

    data["val"][:, :] = spac.value[:, :]
    data.close()
    

def inject_dem():

    args = parse_input_args()
    configurations = load_configuration(args.config)

#-- remap a DEM on to the vertices of the mesh

    print("*inject-dem...")

        # Load the DEM file from the config
    dem_file = configurations['dem_file']

    data = nc.Dataset(dem_file,"r")

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

    """
    # delete groups of dry if too small
    nprt, part = connected_components(
        conn, directed=False, return_labels=True)

    tris = np.argwhere(filt).ravel()
    for iprt in range(nprt):
        itri = np.argwhere(part == iprt)
        if (itri.size <= 2): mask[tris[itri]] = True
    """

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

    args = parse_input_args()
    configurations = load_configuration(args.config)

#-- use the remapped elev. to keep ocean cells

    print("*filter-ocn...")

    elev =(mesh.value[mesh.tria3["index"][:, 0]]
         + mesh.value[mesh.tria3["index"][:, 1]]
         + mesh.value[mesh.tria3["index"][:, 2]]
         + mesh.vmids) / 4.0
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
    
    black_sea = configurations['black_sea']
        # Activate regions based on black_sea option
    if black_sea == 1:  # Caspian and Black Sea
        surf[caspian_region] = -9999.0
        surf[blacksea_region] = -9999.0
    elif black_sea == 2:  # Only Caspian Sea
        surf[caspian_region] = -9999.0
    elif black_sea == 3:  # All except Black Sea
        surf[caspian_region] = -9999.0
        # Additional regions
        additional_region = np.logical_and.reduce((
            mesh.smids[:, 1] >= additional_lat_min,
            mesh.smids[:, 1] <= additional_lat_max,
            mesh.smids[:, 0] >= additional_lon_min,
            mesh.smids[:, 0] <= additional_lon_max
        ))
        surf[additional_region] = 200.

        additional_region2 = np.logical_and.reduce((
            mesh.smids[:, 1] >= additional_lat_min2,
            mesh.smids[:, 1] <= additional_lat_max2,
            mesh.smids[:, 0] >= additional_lon_min2,
            mesh.smids[:, 0] <= additional_lon_max2
        ))
        surf[additional_region2] = 300.

        additional_region3 = np.logical_and.reduce((
            mesh.smids[:, 1] >= additional_lat_min3,
            mesh.smids[:, 1] <= additional_lat_max3,
            mesh.smids[:, 0] >= additional_lon_min3,
            mesh.smids[:, 0] <= additional_lon_max3
        ))
        surf[additional_region3] = 400.
        
    keep = elev <= surf  # only keep tri with wet elev
   
    # iterate on dry cells until none "isolated" 
    knum = np.count_nonzero(keep)
    while (True):   
        keep = filter_dry(mesh, keep)
        if (np.count_nonzero(keep) == knum): break
        knum = np.count_nonzero(keep)
    
    mesh.tria3 = mesh.tria3[keep]

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



if (__name__ == "__main__"):

    args = parse_input_args()
    configurations = load_configuration(args.config)

    create_msh()
    inject_dem()
    filter_ocn()
    
    # viz. in eg. paraview
    jigsawpy.savevtk("test.vtk", mesh)
    
    point = mesh.point["coord"]
    point = jigsawpy.R3toS2(geom.radii, point)  # to [lon,lat] in deg
    point*= 180. / np.pi
    depth = np.reshape(-1*mesh.value, (mesh.value.size, 1))
    depth[depth <= 0] = 2
    point = np.hstack((point, depth))  # append elev. as 3rd coord.
    cells = [("triangle", mesh.tria3["index"])]
    tri_data=cells[0][1]+1
    ww3_mesh_file = configurations['ww3_mesh_file']
    write_gmsh_mesh(ww3_mesh_file, point, tri_data)
