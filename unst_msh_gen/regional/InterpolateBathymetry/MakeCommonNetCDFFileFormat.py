
# This is a script to put all gridded bathymetry datasets into the same format.
# The output files should follow the same format as the 2023 version of the 
# Coastal Relief Model (CRM) available at
# https://www.ncei.noaa.gov/products/coastal-relief-model
# Makes files with same variable names, latitude/longitude conventions, and masking.
#
# To run the full script you will need the following files: 
#
# crm_vol6.nc
# crm_vol8.nc
# crm_southak.nc
# pibhmc_bathy_60m_guam.nc
# ngdc_bathy_90m_amsamoa.nc
# hurl_bathy_60m_nwhi.nc
# RTopo_2_0_4_GEBCO_v2023_60sec_pixel.nc
#
# Sources for these files are within the script in comments.

import numpy as np
import netCDF4 as nc
import math
import subprocess

# Example 1: Basic command execution
# This will run 'ls -l' and print the output to the console.
subprocess.run(["ls", "-l"])

#[0->7] variables(dimensions): |S1 crs(), float64 lat(lat), float64 lon(lon), float32 z(lat, lon)
#[8->9]variables(dimensions): float64 x(x), float64 y(y), float32 z(y, x)
#[10]variables(dimensions): float64 lon(num_lon), float64 lat(num_lat), int16 bed_elevation(num_row, num_col)

Zmin=-12000 # deepest valid value (m) - used to backup masking errors

# from https://www.ncei.noaa.gov/products/coastal-relief-model
flTemplate="crm_vol1_2023.nc"
dataT = nc.Dataset(flTemplate,"r")
print(dataT)


"""
# from https://pae-paha.pacioos.hawaii.edu/erddap/info/srtm30plus_v11_bathy/index.html
flin="srtm30plus_v11_bathy.nc"
flout="srtm30plus_v11_bathy.CRMformat.nc"
data0 = nc.Dataset(flin,"r")
x=np.asarray(data0["lon"][:])
y=np.asarray(data0["lat"][:])
z=np.asarray(data0["elev"][:,:])
#z=float(z)
j=np.where(np.isnan(z))
z[j]=-99999
j=np.where(z<Zmin) #Deeper than M. trench
z[j]=-99999

nx=len(x)
ny=len(y)
with nc.Dataset(flout, 'w', format='NETCDF4') as ncout:
        # Create dimensions
    ncout.createDimension('lon', nx)  # Unlimited dimension
    ncout.createDimension('lat', ny)
    
    lon_var=ncout.createVariable('lon', 'f8', ('lon',))
    lon_var.units         = 'degree_east'
    lon_var.long_name     = 'longitude'
    lon_var.standard_name = 'longitude'
    lon_var.axis          = 'lon'
    lon_var[:]=x[:]

    lat_var=ncout.createVariable('lat', 'f8', ('lat',))
    lat_var.units         = 'degree_north'
    lat_var.long_name     = 'latitude'
    lat_var.standard_name = 'latitude'
    lat_var.axis          = 'lat'
    lat_var[:]=y[:]

    z_var=ncout.createVariable('z', 'f4', ('lat','lon'),fill_value    = -99999)
#    z_var._FillValue    = -99999
    z_var.grid_mapping  = 'crs'
    z_var.long_name     = 'z'
    z_var.units         = 'meters'
    z_var.positive      = 'up'
    z_var.standard_name = 'height'
    z_var.vert_crs_name = 'EGM2008'
    z_var. vert_crs_epsg = 'EPSG:3855'
    z_var[:,:]=z[:,:]
    
    ncout.close
"""


# from https://www.ncei.noaa.gov/products/coastal-relief-model
flin="crm_vol6.nc"
flout="crm_vol6_2023.nc"

data0 = nc.Dataset(flin,"r")
x=np.asarray(data0["x"][:])
y=np.asarray(data0["y"][:])
z=np.asarray(data0["z"][:,:])

j=np.where(np.isnan(z))
z[j]=-99999
j=np.where(z<Zmin) #Deeper than M. trench
z[j]=-99999


nx=len(x)
ny=len(y)

with nc.Dataset(flout, 'w', format='NETCDF4') as ncout:
        # Create dimensions
    ncout.createDimension('lon', nx)  # Unlimited dimension
    ncout.createDimension('lat', ny)
    
    lon_var=ncout.createVariable('lon', 'f8', ('lon',))
    lon_var.units         = 'degree_east'
    lon_var.long_name     = 'longitude'
    lon_var.standard_name = 'longitude'
    lon_var.axis          = 'lon'
    lon_var[:]=x[:]

    lat_var=ncout.createVariable('lat', 'f8', ('lat',))
    lat_var.units         = 'degree_north'
    lat_var.long_name     = 'latitude'
    lat_var.standard_name = 'latitude'
    lat_var.axis          = 'lat'
    lat_var[:]=y[:]

    z_var=ncout.createVariable('z', 'f4', ('lat','lon'),fill_value    = -99999)
#    z_var._FillValue    = -99999
    z_var.grid_mapping  = 'crs'
    z_var.long_name     = 'z'
    z_var.units         = 'meters'
    z_var.positive      = 'up'
    z_var.standard_name = 'height'
    z_var.vert_crs_name = 'EGM2008'
    z_var. vert_crs_epsg = 'EPSG:3855'
    z_var[:,:]=z[:,:]
    
    ncout.close


# from https://www.ncei.noaa.gov/products/coastal-relief-model
flin="crm_vol8.nc"
flout="crm_vol8_2023.nc"

data0 = nc.Dataset(flin,"r")
x=np.asarray(data0["x"][:])
y=np.asarray(data0["y"][:])
z=np.asarray(data0["z"][:,:])

j=np.where(np.isnan(z))
z[j]=-99999
j=np.where(z<Zmin) #Deeper than M. trench
z[j]=-99999


nx=len(x)
ny=len(y)
with nc.Dataset(flout, 'w', format='NETCDF4') as ncout:
        # Create dimensions
    ncout.createDimension('lon', nx)  # Unlimited dimension
    ncout.createDimension('lat', ny)
    
    lon_var=ncout.createVariable('lon', 'f8', ('lon',))
    lon_var.units         = 'degree_east'
    lon_var.long_name     = 'longitude'
    lon_var.standard_name = 'longitude'
    lon_var.axis          = 'lon'
    lon_var[:]=x[:]

    lat_var=ncout.createVariable('lat', 'f8', ('lat',))
    lat_var.units         = 'degree_north'
    lat_var.long_name     = 'latitude'
    lat_var.standard_name = 'latitude'
    lat_var.axis          = 'lat'
    lat_var[:]=y[:]

    z_var=ncout.createVariable('z', 'f4', ('lat','lon'),fill_value    = -99999)
#    z_var._FillValue    = -99999
    z_var.grid_mapping  = 'crs'
    z_var.long_name     = 'z'
    z_var.units         = 'meters'
    z_var.positive      = 'up'
    z_var.standard_name = 'height'
    z_var.vert_crs_name = 'EGM2008'
    z_var. vert_crs_epsg = 'EPSG:3855'
    z_var[:,:]=z[:,:]
    
    ncout.close
    
# from https://pae-paha.pacioos.hawaii.edu/erddap/info/pibhmc_bathy_60m_guam/index.html
flin="pibhmc_bathy_60m_guam.nc"
flout="pibhmc_bathy_60m_guam.crm.nc"
data0 = nc.Dataset(flin,"r")
x=np.asarray(data0["lon"][:])
y=np.asarray(data0["lat"][:])
z=np.asarray(data0["elev"][:,:])

j=np.where(np.isnan(z))
z[j]=-99999
j=np.where(z<Zmin) #Deeper than M. trench
z[j]=-99999

nx=len(x)
ny=len(y)
x=x-360.
with nc.Dataset(flout, 'w', format='NETCDF4') as ncout:
        # Create dimensions
    ncout.createDimension('lon', nx)  # Unlimited dimension
    ncout.createDimension('lat', ny)
    
    lon_var=ncout.createVariable('lon', 'f8', ('lon',))
    lon_var.units         = 'degree_east'
    lon_var.long_name     = 'longitude'
    lon_var.standard_name = 'longitude'
    lon_var.axis          = 'lon'
    lon_var[:]=x[:]

    lat_var=ncout.createVariable('lat', 'f8', ('lat',))
    lat_var.units         = 'degree_north'
    lat_var.long_name     = 'latitude'
    lat_var.standard_name = 'latitude'
    lat_var.axis          = 'lat'
    lat_var[:]=y[:]

    z_var=ncout.createVariable('z', 'f4', ('lat','lon'),fill_value    = -99999)
#    z_var._FillValue    = -99999
    z_var.grid_mapping  = 'crs'
    z_var.long_name     = 'z'
    z_var.units         = 'meters'
    z_var.positive      = 'up'
    z_var.standard_name = 'height'
    z_var.vert_crs_name = 'EGM2008'
    z_var. vert_crs_epsg = 'EPSG:3855'
    z_var[:,:]=z[:,:]
    
    ncout.close
 
# from https://pae-paha.pacioos.hawaii.edu/thredds/bathymetry.html?dataset=ngdc_bathy_90m_amsamoa
flin="ngdc_bathy_90m_amsamoa.nc"
flout="ngdc_bathy_90m_amsamoa.crm.nc"
data0 = nc.Dataset(flin,"r")
x=np.asarray(data0["lon"][:])
y=np.asarray(data0["lat"][:])
z=np.asarray(data0["elev"][:,:])

j=np.where(np.isnan(z))
z[j]=-99999
j=np.where(z<Zmin) #Deeper than M. trench
z[j]=-99999

nx=len(x)
ny=len(y)
with nc.Dataset(flout, 'w', format='NETCDF4') as ncout:
        # Create dimensions
    ncout.createDimension('lon', nx)  # Unlimited dimension
    ncout.createDimension('lat', ny)
    
    lon_var=ncout.createVariable('lon', 'f8', ('lon',))
    lon_var.units         = 'degree_east'
    lon_var.long_name     = 'longitude'
    lon_var.standard_name = 'longitude'
    lon_var.axis          = 'lon'
    lon_var[:]=x[:]

    lat_var=ncout.createVariable('lat', 'f8', ('lat',))
    lat_var.units         = 'degree_north'
    lat_var.long_name     = 'latitude'
    lat_var.standard_name = 'latitude'
    lat_var.axis          = 'lat'
    lat_var[:]=y[:]

    z_var=ncout.createVariable('z', 'f4', ('lat','lon'),fill_value    = -99999)
#    z_var._FillValue    = -99999
    z_var.grid_mapping  = 'crs'
    z_var.long_name     = 'z'
    z_var.units         = 'meters'
    z_var.positive      = 'up'
    z_var.standard_name = 'height'
    z_var.vert_crs_name = 'EGM2008'
    z_var. vert_crs_epsg = 'EPSG:3855'
    z_var[:,:]=z[:,:]
    
    ncout.close

# from https://www.ncei.noaa.gov/metadata/geoportal/rest/metadata/item/gov.noaa.ngdc.mgg.dem:703/html
flin="crm_southak.nc"
flout="crm_southak.CRMformat.nc"
data0 = nc.Dataset(flin,"r")
x=np.asarray(data0["lon"][:])
y=np.asarray(data0["lat"][:])
z=np.asarray(data0["z"][:,:])

j=np.where(z==-2147483648)#Stated fill value
z[j]=-99999
j=np.where(np.isnan(z))
z[j]=-99999
j=np.where(z<Zmin) #Deeper than M. trench
z[j]=-99999

x=x-360
nx=len(x)
ny=len(y)
with nc.Dataset(flout, 'w', format='NETCDF4') as ncout:
        # Create dimensions
    ncout.createDimension('lon', nx)  # Unlimited dimension
    ncout.createDimension('lat', ny)
    
    lon_var=ncout.createVariable('lon', 'f8', ('lon',))
    lon_var.units         = 'degree_east'
    lon_var.long_name     = 'longitude'
    lon_var.standard_name = 'longitude'
    lon_var.axis          = 'lon'
    lon_var[:]=x[:]

    lat_var=ncout.createVariable('lat', 'f8', ('lat',))
    lat_var.units         = 'degree_north'
    lat_var.long_name     = 'latitude'
    lat_var.standard_name = 'latitude'
    lat_var.axis          = 'lat'
    lat_var[:]=y[:]

    z_var=ncout.createVariable('z', 'f4', ('lat','lon'),fill_value    = -99999)
#    z_var._FillValue    = -99999
    z_var.grid_mapping  = 'crs'
    z_var.long_name     = 'z'
    z_var.units         = 'meters'
    z_var.positive      = 'up'
    z_var.standard_name = 'height'
    z_var.vert_crs_name = 'EGM2008'
    z_var. vert_crs_epsg = 'EPSG:3855'
    z_var[:,:]=z[:,:]
    
    ncout.close

# from https://pae-paha.pacioos.hawaii.edu/thredds/bathymetry.html?dataset=hurl_bathy_60m_nwhi
flin="hurl_bathy_60m_nwhi.nc"
flout="hurl_bathy_60m_nwhi.CRMformat.nc"
data0 = nc.Dataset(flin,"r")
x=np.asarray(data0["lon"][:])
y=np.asarray(data0["lat"][:])
z=np.asarray(data0["elev"][:,:])
#z=float(z)
j=np.where(z==-2147483648)#Stated fill value
z[j]=-99999
j=np.where(np.isnan(z))
z[j]=-99999
j=np.where(z<Zmin) #Deeper than M. trench
z[j]=-99999
#x=x-360
nx=len(x)
ny=len(y)
with nc.Dataset(flout, 'w', format='NETCDF4') as ncout:
        # Create dimensions
    ncout.createDimension('lon', nx)  # Unlimited dimension
    ncout.createDimension('lat', ny)
    
    lon_var=ncout.createVariable('lon', 'f8', ('lon',))
    lon_var.units         = 'degree_east'
    lon_var.long_name     = 'longitude'
    lon_var.standard_name = 'longitude'
    lon_var.axis          = 'lon'
    lon_var[:]=x[:]

    lat_var=ncout.createVariable('lat', 'f8', ('lat',))
    lat_var.units         = 'degree_north'
    lat_var.long_name     = 'latitude'
    lat_var.standard_name = 'latitude'
    lat_var.axis          = 'lat'
    lat_var[:]=y[:]

    z_var=ncout.createVariable('z', 'f4', ('lat','lon'),fill_value    = -99999)
#    z_var._FillValue    = -99999
    z_var.grid_mapping  = 'crs'
    z_var.long_name     = 'z'
    z_var.units         = 'meters'
    z_var.positive      = 'up'
    z_var.standard_name = 'height'
    z_var.vert_crs_name = 'EGM2008'
    z_var. vert_crs_epsg = 'EPSG:3855'
    z_var[:,:]=z[:,:]
   
    ncout.close

#Global gridded bathymetry data
# $wget https://github.com/dengwirda/dem/releases/download/v0.1.1/RTopo_2_0_4_GEBCO_v2023_60sec_pixel.zip

flin="RTopo_2_0_4_GEBCO_v2023_60sec_pixel.nc"
flout="RTopo_2_0_4_GEBCO_v2023_60sec_pixel.CRMformat.nc"

data0 = nc.Dataset(flin,"r")
x0=np.asarray(data0["lon"][:])
y0=np.asarray(data0["lat"][:])
z=np.asarray(data0["bed_elevation"][:,:])
nx0=len(x0)
ny0=len(y0)
nx=nx0-1
ny=ny0-1
x=np.zeros(nx)
y=np.zeros(ny)
for k in range(nx):
    x[k]=(x0[k]+x0[k+1])/2
for k in range(ny):
    y[k]=(y0[k]+y0[k+1])/2
    
with nc.Dataset(flout, 'w', format='NETCDF4') as ncout:
        # Create dimensions
    ncout.createDimension('lon', nx)  # Unlimited dimension
    ncout.createDimension('lat', ny)
    
    lon_var=ncout.createVariable('lon', 'f8', ('lon',))
    lon_var.units         = 'degree_east'
    lon_var.long_name     = 'longitude'
    lon_var.standard_name = 'longitude'
    lon_var.axis          = 'lon'
    lon_var[:]=x[:]

    lat_var=ncout.createVariable('lat', 'f8', ('lat',))
    lat_var.units         = 'degree_north'
    lat_var.long_name     = 'latitude'
    lat_var.standard_name = 'latitude'
    lat_var.axis          = 'lat'
    lat_var[:]=y[:]

    z_var=ncout.createVariable('z', 'f4', ('lat','lon'),fill_value    = -99999)
    z_var.grid_mapping  = 'crs'
    z_var.long_name     = 'z'
    z_var.units         = 'meters'
    z_var.positive      = 'up'
    z_var.standard_name = 'height'
    z_var.vert_crs_name = 'EGM2008'
    z_var. vert_crs_epsg = 'EPSG:3855'
    z_var[:,:]=z[:,:]
    
    ncout.close

    
