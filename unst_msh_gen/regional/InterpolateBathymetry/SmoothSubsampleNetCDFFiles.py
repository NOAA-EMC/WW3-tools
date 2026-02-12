
# This is a script should be run after MakeCommonNetCDFFileFormat.py
#
# The purpose is to smooth the various high resolution local gridded 
# bathymetry to a usefull length scale for mesh generation(this will
# vary with mesh generation application). The current choices are for
# meshes with a minimum lengthscale of ~500m.
#
# Most of the files are Coastal Relief Model (CRM) 
# https://www.ncei.noaa.gov/products/coastal-relief-model
#
# Other files are mostly derived from from bathymetry available at:
# https://pae-paha.pacioos.hawaii.edu/thredds/bathymetry.html



import numpy as np
import netCDF4 as nc
from scipy import signal

#[0->7] variables(dimensions): |S1 crs(), float64 lat(lat), float64 lon(lon), float32 z(lat, lon)
#[8->9]variables(dimensions): float64 x(x), float64 y(y), float32 z(y, x)
#[10]variables(dimensions): float64 lon(num_lon), float64 lat(num_lat), int16 bed_elevation(num_row, num_col)

flnms=[
    "crm_vol1_2023.nc",#dx= 0.00027, 3600/deg
    "crm_vol2_2023.nc",
    "crm_vol3_2023.nc",
    "crm_vol4_2023.nc",
    "crm_vol5_2023.nc",
    "crm_vol7_2024.nc",
    "crm_vol9_2023.nc",
    "crm_vol10_2023.nc", #problem with mask values
    "crm_vol6_2023.nc",#dx= 0.0008 cut to 3 x 3
    "crm_vol8_2023.nc",#dx= 0.0008,1200 / deg
    "hurl_bathy_60m_nwhi.CRMformat.nc",
    "pibhmc_bathy_60m_guam.crm.nc",
    "ngdc_bathy_90m_amsamoa.crm.nc"
    ]

count=0
filval=-99999
for fl in flnms:
    count=count+1
    
    Nsubsample=9
    Nstart=4
    if (fl=="crm_vol6_2023.nc"or fl=="crm_vol8_2023.nc"  ):
        Nsubsample=3
        Nstart=1
    if fl == "hurl_bathy_60m_nwhi.CRMformat.nc":
        Nsubsample=5
        Nstart=2
    if (fl == "pibhmc_bathy_60m_guam.crm.nc" or "ngdc_bathy_90m_amsamoa.crm.nc"):
        Nsubsample=3
        Nstart=1
        
    data0 = nc.Dataset(fl,"r")
    x=np.asarray(data0["lon"][:])
    y=np.asarray(data0["lat"][:])
    if fl == "crm_socal_1as_vers2.nc":
        z=np.asarray(data0["Band1"][:])
        z=z.astype(np.float32)
    else:
        z=np.asarray(data0["z"][:])
    nx=len(x)
    ny=len(y)
    dx=x[1]-x[0]
    dy=y[1]-y[0]
    print(fl)
    print(z.shape)
    print("org. dx= "+str(dx)) 
    print(1.0/dx) 
    print("org. dy= "+str(dy)) 
    print(1.0/dy)
    if ( (fl == "crm_vol10_2023.nc") or (fl == "crm_socal_1as_vers2.nc") ) : # This file has an un recognized fill value
        print("fixing bad fill value in: "+fl)
        jb=np.where(z==-9999)
        z[jb]= -float("inf")
    
    W=np.ones((Nsubsample,Nsubsample))/(Nsubsample**2)
    z = signal.convolve2d(z,W, mode='same')
    x=x[Nstart :: Nsubsample]
    y=y[Nstart :: Nsubsample]
    z=z[Nstart :: Nsubsample,Nstart :: Nsubsample]
    nx=len(x)
    ny=len(y)
#    flout=fl+'.S250m.nc'
    flout=fl+'.S250m.VB.nc'
        # This file has an un recognized fill value
#    if fl == "crm_vol10_2023.nc": 
    if ( (fl == "crm_vol10_2023.nc") or (fl == "crm_socal_1as_vers2.nc") ) : # This file has an un recognized fill value
        print("patching correct fill values in: "+flout)
        jb=np.where(z<filval)
        z[jb] = filval
    
    with nc.Dataset(flout, 'w', format='NETCDF4') as ncout:
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

        z_var=ncout.createVariable('z', 'f4', ('lat','lon'),fill_value    = filval)
        z_var.grid_mapping  = 'crs'
        z_var.long_name     = 'z'
        z_var.units         = 'meters'
        z_var.positive      = 'up'
        z_var.standard_name = 'height'
        z_var.vert_crs_name = 'EGM2008'
        z_var. vert_crs_epsg = 'EPSG:3855'
        z_var[:,:]=z[:,:]
        
        ncout.close
        
    dx=x[1]-x[0]
    dy=y[1]-y[0]
    print("new dx= "+str(dx)) 
    print(1.0/dx) 
    print("new dy= "+str(dy)) 
    print(1.0/dy) 
    print(z.shape)
    print("--------------------------------------")

