#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
organizeDistanceToCoast.py

VERSION AND LAST UPDATE:
 v1.0  04/04/2022

PURPOSE:
 Small code to organize the GSFC/NASA global information of distance to 
  the coast, into a gridded netcdf file.

USAGE:
 Download dist2coast.txt from
  https://oceancolor.gsfc.nasa.gov/docs/distfromcoast/
 And define lat and lon arrays where you want the information,
  rlon, rlon (edit below)

OUTPUT:
 netcdf file distFromCoast.nc with distance to the nearest coast
  at rlon rlat points defined by the user.

DEPENDENCIES:
 See dependencies.py and the imports below.

AUTHOR and DATE:
 04/04/2022: Ricardo M. Campos, first version.

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import numpy as np
from matplotlib.mlab import *
from pylab import *
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import cm
from netCDF4 import Dataset
import warnings; warnings.filterwarnings("ignore")
colormap = cm.GMT_polar
palette = plt.cm.jet
palette.set_bad('aqua', 10.0)
# netcdf format
fnetcdf="NETCDF4"

# GSFC/NASA file can be downloaded at
# https://oceancolor.gsfc.nasa.gov/docs/distfromcoast/
rdfc=np.loadtxt('dist2coast.txt', usecols=(2,))
rlon=np.arange(-179.98,179.98,0.04)
rlat=np.arange(-89.98,89.98,0.04)

dfc = np.zeros((rlat.shape[0],rlon.shape[0]),'f')

c=0
for i in xrange(0,rlat.shape[0]):
	for j in xrange(0,rlon.shape[0]):
		dfc[-(i+1),j] = rdfc[c]
		c = c+1
	print(repr(i))

# open a new netCDF file for writing.
ncfile = Dataset('distFromCoast.nc', "w", format=fnetcdf) 
ncfile.history="Distance to coast GSFC/NASA and NOAA" 
# create the lat and lon dimensions.
ncfile.createDimension( 'latitude' , rlat.shape[0] ) 
ncfile.createDimension( 'longitude' , rlon.shape[0] ) 
lats = ncfile.createVariable('latitude',dtype('float32').char,('latitude',)) 
lons = ncfile.createVariable('longitude',dtype('float32').char,('longitude',)) 
# Assign units attributes to coordinate var data. This attaches a text attribute to each of the coordinate variables, containing the units.
lats.units = 'degrees_north'
lons.units = 'degrees_east'
fdfc.units = 'km'
# write data to coordinate vars.
lats[:] = rlat
lons[:] = rlon
# create  variable
fdfc = ncfile.createVariable('distcoast',dtype('float32').char,('latitude','longitude'))
# write data to variables.
fdfc[:,:]=dfc[:,:] 
# close the file
ncfile.close()
print('netcdf ')

