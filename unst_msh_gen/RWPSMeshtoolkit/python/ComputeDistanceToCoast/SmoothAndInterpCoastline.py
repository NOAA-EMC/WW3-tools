import geopandas as gpd
import time
from SmoothAndSubsampleCoastline import *

import numpy as np
import netCDF4 as nc
import jigsawpy
import geopandas as gpd


fl='../RWPS/Data/GlobalCoast/GSHHS_shp/f/GSHHS_f_L1.shp'
    #gdf = gpd.read_file("../RWPS/Data/us_coastline/tl_2023_us_coastline.shp")
dxS=10000.       #(positive real)Smooth coastline to dxS meters length scale 
dxI=5000.       #(positive real)Interpolate smoothed coastline to dxI meters
edges=SmoothAndSubsampleCoastlineGeom(fl,dxS,dxI,100)
lonUS=points[:,0]
latUS=points[:,1]
