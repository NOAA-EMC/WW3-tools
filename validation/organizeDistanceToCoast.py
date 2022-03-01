import numpy as np
from matplotlib.mlab import *
from pylab import *
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import cm
from netCDF4 import Dataset
colormap = cm.GMT_polar
palette = plt.cm.jet
palette.set_bad('aqua', 10.0)

# GSFC/NASA file can be downloaded at
# https://oceancolor.gsfc.nasa.gov/docs/distfromcoast/
rdfc=np.loadtxt('dist2coast.txt', usecols=(2,))
rlon=frange(-179.98,179.98,0.04)
rlat=frange(-89.98,89.98,0.04)

dfc = np.zeros((rlat.shape[0],rlon.shape[0]),'f')

c=0
for i in xrange(0,rlat.shape[0]):
	for j in xrange(0,rlon.shape[0]):
		dfc[-(i+1),j] = rdfc[c]
		c = c+1
	print(repr(i))

#save netcdf
fnetcdf="NETCDF3_CLASSIC"
# open a new netCDF file for writing.
ncfile = Dataset('distFromCoast.nc', "w", format=fnetcdf) 
ncfile.history="Distance to coast NASA/NOAA" 
# create the lat and lon dimensions.
ncfile.createDimension( 'latitude' , rlat.shape[0] ) 
ncfile.createDimension( 'longitude' , rlon.shape[0] ) 
lats = ncfile.createVariable('latitude',dtype('float32').char,('latitude',)) 
lons = ncfile.createVariable('longitude',dtype('float32').char,('longitude',)) 
# Assign units attributes to coordinate var data. This attaches a text attribute to each of the coordinate variables, containing the units.
lats.units = 'degrees_north'
lons.units = 'degrees_east'
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
