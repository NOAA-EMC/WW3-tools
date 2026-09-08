from datetime import datetime

import bufr
import numpy as np
import netCDF4 as nc
import pandas as pd

def main():
    '''
    Convert sofar data from bufr to a netcdf file.
    Reference for bufr data tables: https://www.nco.ncep.noaa.gov/sib/jeff/bufrtab_tableb.html
    '''
    # Set up the file path to the 
    # file_path = './data/bufr_data/sofar_b001_xx121_20260101'
    file_path = './data/bufr_data/sofar_b001_xx121_20260102'
    nc_path = './data/nc_data/sofar_b001_xx121_20260102.nc' # ADJUST to loop

    # Set up the variables to query
    q = bufr.QuerySet()
    #q.add('report_identifier', '*/RPID')

    # Time Variables
    q.add('year', '*/YEAR')
    q.add('month', '*/MNTH')
    q.add('day', '*/DAYS')
    q.add('hour', '*/HOUR')
    q.add('minute', '*/MINU')

    # Identifier Variables
    q.add('station_name', '*/LSTN')
    # q.add('buoy_id', '*/BPID') - Empty
    q.add('buoy_type', '*/BUYT')

    # Location Variables
    q.add('latitude', '*/CLATH')
    q.add('longitude', '*/CLONH')
   

    # Data Variables
    q.add('hs', '*/BASICWAV/SGWH') # Significant Wave Height
    q.add('dsdw', '*/BASICWAV/DSDW') # Directional Spread of dominant wave
    q.add('ddwc', '*/BASICWAV/DDWC') # Direction from which dominant waves are coming
    q.add('avwp', '*/BASICWAV/AVWP') # Average Wave Period
    q.add('spwp', '*/BASICWAV/SPWP') # Spectral Peak Wave Period
    q.add('wspd', '*/WSPD') # Wind Speed
    q.add('wdir', '*/WDIR') # Wind Direction
    # q.add('duration_of_record', '*/BASICWAV/DOWR') - Empty
    # q.add('h_max', '*/BASICWAV/MXWH') - Empty
    
    # Open the bufr file
    with bufr.File(file_path) as f:
        r = f.execute(q)

    # Open the netcdf file for writing 
    rootgrp = nc.Dataset(nc_path, 'w', clobber=True)

    # netCDF metadata
    rootgrp.description = "Example sofar buoy data bufr"
    rootgrp.source = "Converted from NCO BUFR data tanks"

    # Time
    year = r.get('year')
    month = r.get('month')
    day = r.get('day')
    hour = r.get('hour')
    minute = r.get('minute')
    time_dim = rootgrp.createDimension('time', year.size)
    time_nc = rootgrp.createVariable('time', 'f8', ('time',))
    time_nc.units = "seconds since 1970-01-01 00:00:00"
    time_nc.calendar = "gregorian"

    # Create datetime, convert to epoch time and save to nc
    datetime_array = np.empty(shape=year.shape)
    for n in np.arange(year.size): 
        dt = datetime(year[n], month[n], day[n], hour[n], minute[n])
        datetime_array[n] = nc.date2num(dt, units=time_nc.units,calendar=time_nc.calendar)

    # Save the epoch time to the netCDF
    time_nc[:] = datetime_array

    # Latitude
    lats = rootgrp.createVariable('lat', 'f4', ('time',))
    lats.long_name = 'latitude'
    lats.units = 'degrees_north'
    lats[:] = r.get('latitude')

    # Longitude
    lons = rootgrp.createVariable('lon', 'f4', ('time',))
    lons.long_name = 'longitude'
    lons.units = 'degrees_east'
    lons[:] = r.get('longitude')

    # Long Station name - Buoy Identifier
    lstn = rootgrp.createVariable('lstn', str, ('time',))
    lstn.long_name = 'long station name'
    lstn.units = 'str'
    lstn[:] = r.get('station_name')

    # Buoy Type
    buoy_type = rootgrp.createVariable('buoy_type', 'i4', ('time',))
    buoy_type.units = 'int'
    buoy_type.comment = 'Type 9 is a SOFAR Buoy, ref: https://www.nco.ncep.noaa.gov/sib/jeff/CodeFlag_0_STDv31_LOC7.html'
    buoy_type[:] = r.get('buoy_type')

    # Significant Wave Height
    hs = rootgrp.createVariable('hs', 'f4', ('time',))
    hs.long_name = 'significant wave height'
    hs.units = 'meters'
    hs[:] = r.get('hs')

    # Directional Spread of Dominant Wave
    dsdw = rootgrp.createVariable('dsdw', 'i4', ('time',))
    dsdw.long_name = 'directional spread of dominant wave'
    dsdw.units = 'degree'
    dsdw[:] = r.get('dsdw')

    # Direction from which dominant waves are coming
    ddwc = rootgrp.createVariable('ddwc', 'i4', ('time',))
    ddwc.long_name = 'direction from which dominant waves are coming'
    ddwc.units = 'degree true'
    ddwc[:] = r.get('ddwc')

    # Average Wave Period
    avwp = rootgrp.createVariable('avwp', 'f4', ('time',))
    avwp.long_name = 'average wave period'
    avwp.units = 'seconds'
    avwp[:] = r.get('avwp')

    # Spectral peak wave period
    spwp = rootgrp.createVariable('spwp', 'f4', ('time',))
    spwp.long_name = 'spectral peak wave period'
    spwp.units = 'seconds'
    spwp[:] = r.get('spwp')

    # Wind Speed
    wspd = rootgrp.createVariable('wspd', 'f4', ('time',))
    wspd.long_name = 'wind speed'
    wspd.units = 'meters per second'
    wspd[:] = r.get('wspd')

    # Wind Direction
    wdir = rootgrp.createVariable('wdir', 'f4', ('time',))
    wdir.long_name = 'wind direction'
    wdir.units = 'degree true'
    wdir[:] = r.get('wdir')

    # Maximum Wave Height - Empty
    # h_max = rootgrp.createVariable('h_max', 'f4', ('time',))
    # h_max.units = 'meters'
    # h_max[:] = r.get('h_max')
    # print(r.get('h_max'))
    # # 

        
    #     # # 5. Add variable-specific metadata (CF Conventions)
    #     # temp = nc_file.createVariable('temperature', 'f4', ('time', 'lat', 'lon',))
    #     # temp.units = 'Kelvin'
    #     # temp.long_name = 'Surface Air Temperature'
    #     # temp_data = np.random.uniform(low=270, high=300, size=(3, 4, 5))
    #     # temp[:, :, :] = temp_data   


    return

if __name__=="__main__":
    main()