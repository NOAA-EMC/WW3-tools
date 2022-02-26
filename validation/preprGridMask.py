import numpy as np
from matplotlib.mlab import *
from pylab import *
import matplotlib
# matplotlib.use('Agg')
import xarray as xr
import netCDF4 as nc
from pylab import *
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
from matplotlib import ticker
# import pickle
import warnings
warnings.filterwarnings("ignore")

palette = plt.cm.jet
# Font size and style
sl=14
matplotlib.rcParams.update({'font.size': sl}); plt.rc('font', size=sl) 
matplotlib.rc('xtick', labelsize=sl); matplotlib.rc('ytick', labelsize=sl); matplotlib.rcParams.update({'font.size': sl})

# -----------------
# Grid Name
gridn='GEFS'
# Minimum water depth
mindepth=80. # in meters
# Minimum distance from the coast
mindfc=50. # in Km
# -----------------

# Sample file Model(reference) lat lon array
ds=xr.open_dataset('ww3gefs.20160921_field.nc')
mapsta=ds["MAPSTA"].values[:,:]
lat = ds['latitude'].values[:]; lon = ds['longitude'].values[:]
ds.close(); del ds

# BATHYMETRY Etopo grid
ds = xr.open_dataset('etopo1.nc')
latb = ds['lat'].values[:]; lonb = ds['lon'].values[:]
b = ds['z'].values[:,:]
# interpolate Bathymetry to Model
dsi = ds.interp(lat=lat[:], lon=lon[:])
ib = dsi['z'].values[:,:]
ds.close(); del ds, dsi

# distance to the coast
ds = xr.open_dataset('distFromCoast.nc')
latd = ds['latitude'].values[:]; lond = ds['longitude'].values[:]
dfc = ds['distcoast'].values[:,:]
# interpolate to Model
dsc = ds.interp(latitude=lat[:], longitude=lon[:])
idfc = dsc['distcoast'].values[:,:]; ds.close(); del ds, dsc

# Build Mask (nan = land excluded; 0 = ocean excluded; 1 = ocean valid)
mask = np.zeros((lat.shape[0],lon.shape[0]),'f')
# excluding continent or GEFS mask
ind = np.where((ib>0)|(mapsta>100)|(mapsta==0)); mask[ind[0],ind[1]] = np.nan; del ind
# excluding depth dist-to-coast
ind = np.where( (ib<=(-1*mindepth)) & (idfc>=mindfc) & (np.isnan(mask)==False) ); mask[ind[0],ind[1]] = 1; del ind


# Plots -----
# Water depth is positive, by definition
ib = np.array(ib*-1); ib[ib<0]=np.nan
#
ib[ib<0.]=np.nan; ib[np.isnan(mask)==True]=np.nan
levels = np.linspace( ib[(np.isnan(ib)==False)].min(), np.percentile(ib[(np.isnan(ib)==False)],99.), 100)
#
fig, ax = plt.subplots(figsize=(7,4))
# ax = plt.axes(projection=ccrs.Robinson())
ax = plt.axes(projection=ccrs.PlateCarree()) 
ax.set_extent([lon.min(),lon.max(),lat.min(),lat.max()], crs=ccrs.PlateCarree())
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
plt.contourf(lon,lat,ib,levels,transform = ccrs.PlateCarree(),cmap=palette,extend="max", zorder=2)
ax.add_feature(cartopy.feature.OCEAN,facecolor=("white"))
ax.add_feature(cartopy.feature.LAND,facecolor=("lightgrey"), edgecolor='grey',linewidth=0.5)
ax.add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1)
ax.coastlines(resolution='110m', color='grey',linewidth=0.5, linestyle='-', alpha=1)
fig.tight_layout()
ax = plt.gca()
pos = ax.get_position()
l, b, w, h = pos.bounds
cax = plt.axes([l+0.07, b-0.07, w-0.15, 0.025]) # setup colorbar axes.
cbar=plt.colorbar(cax=cax, orientation='horizontal'); cbar.ax.tick_params(labelsize=10)
tick_locator = ticker.MaxNLocator(nbins=7); cbar.locator = tick_locator; cbar.update_ticks()
plt.axes(ax)  # make the original axes current again
fig.tight_layout()
# plt.savefig('wfields_'+wvar+'_'+np.str(pd.to_datetime(wtime[::sk][t]).strftime('%Y%m%d%H'))+'.eps', format='eps', dpi=200)
plt.savefig('bathymetry_'+gridn+'.png', dpi=300, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
# pickle.dump(fig, open('wfields_'+wvar+'_'+np.str(pd.to_datetime(wtime[::sk][t]).strftime('%Y%m%d%H'))+'.pickle', 'wb'))
plt.close('all'); del ax, fig


idfc[ib<0.]=NaN; idfc[np.isnan(mask)==True]=np.nan
levels = np.linspace( idfc[(np.isnan(idfc)==False)].min(), np.percentile(idfc[(np.isnan(idfc)==False)],99.), 100)

fig, ax = plt.subplots(figsize=(7,4))
# ax = plt.axes(projection=ccrs.Robinson())
ax = plt.axes(projection=ccrs.PlateCarree()) 
ax.set_extent([lon.min(),lon.max(),lat.min(),lat.max()], crs=ccrs.PlateCarree())
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
plt.contourf(lon,lat,idfc,levels,transform = ccrs.PlateCarree(),cmap=palette,extend="max", zorder=2)
ax.add_feature(cartopy.feature.OCEAN,facecolor=("white"))
ax.add_feature(cartopy.feature.LAND,facecolor=("lightgrey"), edgecolor='grey',linewidth=0.5)
ax.add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1)
ax.coastlines(resolution='110m', color='grey',linewidth=0.5, linestyle='-', alpha=1)
fig.tight_layout()
ax = plt.gca()
pos = ax.get_position()
l, b, w, h = pos.bounds
cax = plt.axes([l+0.07, b-0.07, w-0.15, 0.025]) # setup colorbar axes.
cbar=plt.colorbar(cax=cax, orientation='horizontal'); cbar.ax.tick_params(labelsize=10)
tick_locator = ticker.MaxNLocator(nbins=7); cbar.locator = tick_locator; cbar.update_ticks()
plt.axes(ax)  # make the original axes current again
fig.tight_layout()
# plt.savefig('wfields_'+wvar+'_'+np.str(pd.to_datetime(wtime[::sk][t]).strftime('%Y%m%d%H'))+'.eps', format='eps', dpi=200)
plt.savefig('distanceToCoast_'+gridn+'.png', dpi=300, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
# pickle.dump(fig, open('wfields_'+wvar+'_'+np.str(pd.to_datetime(wtime[::sk][t]).strftime('%Y%m%d%H'))+'.pickle', 'wb'))
plt.close('all'); del ax, fig

mask[:,0]=mask[:,-1]
levels = np.linspace(-2,3,10)
fig, ax = plt.subplots(figsize=(7,4))
# ax = plt.axes(projection=ccrs.Robinson())
ax = plt.axes(projection=ccrs.PlateCarree()) 
ax.set_extent([lon.min(),lon.max(),lat.min(),lat.max()], crs=ccrs.PlateCarree())
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
plt.contourf(lon,lat,-mask,levels,transform = ccrs.PlateCarree(),cmap=palette,extend="max", zorder=2)
ax.add_feature(cartopy.feature.OCEAN,facecolor=("white"), zorder=1)
ax.add_feature(cartopy.feature.LAND,facecolor=("lightgrey"), edgecolor='grey',linewidth=0.5, zorder=3)
ax.add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=2, zorder=3)
ax.coastlines(resolution='110m', color='grey',linewidth=0.5, linestyle='-', alpha=1, zorder=3)
fig.tight_layout()
# plt.savefig('wfields_'+wvar+'_'+np.str(pd.to_datetime(wtime[::sk][t]).strftime('%Y%m%d%H'))+'.eps', format='eps', dpi=200)
plt.savefig('Mask_'+gridn+'.png', dpi=300, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
# pickle.dump(fig, open('wfields_'+wvar+'_'+np.str(pd.to_datetime(wtime[::sk][t]).strftime('%Y%m%d%H'))+'.pickle', 'wb'))
plt.close('all'); del ax, fig

print('Number of Ocean points: '+repr(mask[mask>=0].shape[0]))
print('Number of Ocean points valid: '+repr(mask[mask>0].shape[0]))

#save netcdf ==================
fnetcdf="NETCDF4"
# open a new netCDF file for writing.
ncfile = nc.Dataset('gridInfo_'+gridn+'.nc', "w", format=fnetcdf) 
ncfile.description='Bathymetry, Distance from the coast, and Mask. Total of '+repr(mask[mask>=0].shape[0])+' Ocean grid points, and '+repr(mask[mask>0].shape[0])+' valid ocean grid points to use.'
# create the lat and lon dimensions.
ncfile.createDimension( 'latitude' , lat.shape[0] ) 
ncfile.createDimension( 'longitude' , lon.shape[0] ) 
lats = ncfile.createVariable('latitude',dtype('float32').char,('latitude',)) 
lons = ncfile.createVariable('longitude',dtype('float32').char,('longitude',)) 
# Assign units attributes to coordinate var data. This attaches a text attribute to each of the coordinate variables, containing the units.
lats.units = 'degrees_north'
lons.units = 'degrees_east'
# write data to coordinate vars.
lats[:] = lat
lons[:] = lon
# create  variable
fdfc = ncfile.createVariable('distcoast',dtype('float32').char,('latitude','longitude'))
fbt = ncfile.createVariable('depth',dtype('float32').char,('latitude','longitude'))
fm = ncfile.createVariable('mask',dtype('float32').char,('latitude','longitude'))
# write data to variables.
fdfc[:,:]=idfc[:,:] 
fbt[:,:]=ib[:,:] 
fm[:,:]=mask[:,:]
# close the file
ncfile.close()
print('netcdf ok')

