function [filename] = write_directional_spectra_nc(filename,testcase,...
    pointID,Lat,Lon,dpt,wndspd,wnddir,curspd,curdir,time,freq,dir,EFTH)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function writes Directional Spectral density file in %
% WW3 netcdf format                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Ali Abdolali April 2017 ali.abdolali@noaa.gov        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testcase name: i.e. Inlet test'
% time
% ntime: length of time
% nfreq: number of frequency band
% nDir: number of directional band
% pointnumber : number of points- use 1 for this variable
% freq: frequency (Hz)
% dir: direction (rad)
% pointID: the boundary point name: i.e. 'b42001' 
% Lat: latitude (degrees)
% Lon: longitude (degrees)
% dpt: depth (m) [pointnumber x 1]
% wndspd: wind  speed (m/s) [pointnumber, ntime]
% wnddir: wind direction (degrees) [pointnumber, ntime]
% curspd: current speed (m/s) [pointnumber, ntime]
% curdir: current direction (degrees) [pointnumber, ntime]
% EFTH: directional spectral density [nDir x nfreq  x pointnumber x ntime]
% filename: name of output file 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%example: [filename] = write_directional_spectra_nc('B42001.nc',...
%'inlet','B42001',42,221.1,1521,30,270,0.5,36,time,freq,dir,EFTH)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nDir,nfreq,nStation,ntime] = size(EFTH);
frequency1=freq*1.047619000313776;
frequency1(1)=1;
frequency2=freq*0.952380928728317;
frequency2(end)=1;

nc = netcdf.create(ncfile, '64BIT_OFFSET');

% define global attributes
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'start_date',datestr(time(1)))
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'stop_date',datestr(time(end)))
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'source',testcase)
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'field_type','n/a')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'data_type','OCO spectra 2D')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'southernmost_latitude','n/a')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'northernmost_latitude','n/a')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'latitude_resolution','n/a')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'westernmost_longitude','n/a')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'easternmost_longitude','n/a')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'longitude_resolution','n/a')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'minimum_altitude','n/a')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'maximum_altitude','n/a')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'altitude_resolution','n/a')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'area','Global')
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'format_version','1.1')
                  

% define dimensions

unlim_id = netcdf.defDim(nc, 'time', netcdf.getConstant('NC_UNLIMITED'));
STAT = netcdf.defDim(nc, 'station', nStation);
STRING = netcdf.defDim(nc, 'string16', 16);
FREQ = netcdf.defDim(nc, 'frequency', nfreq);
DIR = netcdf.defDim(nc, 'direction', nDir);

STATION=pointID;

time_varid = netcdf.defVar(nc, 'time', 'double', unlim_id);
netcdf.putAtt(nc, time_varid, 'long_name', 'julian day (UT)');
netcdf.putAtt(nc, time_varid, 'units', 'days since 1990-01-01 00:00:00');
netcdf.putAtt(nc, time_varid, 'field', 'time');
netcdf.putAtt(nc, time_varid, 'conventions', 'Relative julian days with decimal part (as parts of the day )');
netcdf.putAtt(nc, time_varid, 'axis', 'T');

station_varid = netcdf.defVar(nc, 'station', 'NC_INT', [STAT]);
netcdf.putAtt(nc, station_varid, 'long_name', 'station id');
netcdf.putAtt(nc, station_varid, 'axis', 'X');

station16_varid = netcdf.defVar(nc, 'string16', 'NC_INT', [STRING]);
netcdf.putAtt(nc, station_varid, 'long_name', 'station_name number of characters');
netcdf.putAtt(nc, station_varid, 'axis', 'W');

station_name_varid = netcdf.defVar(nc, 'station_name', 'NC_CHAR', [STAT, STRING]);
netcdf.putAtt(nc, station_varid, 'long_name', 'station name');
netcdf.putAtt(nc, station_varid, 'axis', 'XW');

                 
lon_varid=netcdf.defVar(nc, 'x' ,'NC_FLOAT',[STAT,unlim_id]);
netcdf.putAtt(nc, lon_varid, 'long_name', 'x');
netcdf.putAtt(nc, lon_varid, 'standard_name', 'x');
netcdf.putAtt(nc, lon_varid, 'globwave_name', 'x');
netcdf.putAtt(nc, lon_varid, 'content', 'TX');
netcdf.putAtt(nc, lon_varid, 'associates', 'time station');
netcdf.putAtt(nc, lon_varid, 'units', 'm');

lat_varid=netcdf.defVar(nc, 'y' ,'NC_FLOAT',[STAT,unlim_id]);
netcdf.putAtt(nc, lat_varid, 'long_name', 'y');
netcdf.putAtt(nc, lat_varid, 'standard_name', 'y');
netcdf.putAtt(nc, lat_varid, 'globwave_name', 'y');
netcdf.putAtt(nc, lat_varid, 'content', 'TX');
netcdf.putAtt(nc, lat_varid, 'associates', 'time station');
netcdf.putAtt(nc, lat_varid, 'units', 'm');


f_varid=netcdf.defVar(nc, 'frequency' ,'NC_FLOAT',[FREQ]);
netcdf.putAtt(nc, f_varid, 'long_name', 'frequency of center band');
netcdf.putAtt(nc, f_varid, 'standard_name', 'sea_surface_wave_frequency');
netcdf.putAtt(nc, f_varid, 'globwave_name', 'frequency');
netcdf.putAtt(nc, f_varid, 'units', 's-1');
netcdf.putAtt(nc, f_varid, 'axis', 'Y');

f1_varid=netcdf.defVar(nc, 'frequency1' ,'NC_FLOAT',[FREQ]);
netcdf.putAtt(nc, f1_varid, 'long_name', 'frequency of lower band');
netcdf.putAtt(nc, f1_varid, 'standard_name', 'frequency_of_lower_band');
netcdf.putAtt(nc, f1_varid, 'globwave_name', 'frequency_lower_band');
netcdf.putAtt(nc, f1_varid, 'units', 's-1');
netcdf.putAtt(nc, f1_varid, 'axis', 'Y');
netcdf.putAtt(nc, f1_varid, 'associates', 'frequency');

f2_varid=netcdf.defVar(nc, 'frequency2' ,'NC_FLOAT',[FREQ]);
netcdf.putAtt(nc, f2_varid, 'long_name', 'frequency of upper band');
netcdf.putAtt(nc, f2_varid, 'standard_name', 'frequency_of_upper_band');
netcdf.putAtt(nc, f2_varid, 'globwave_name', 'frequency_upper_band');
netcdf.putAtt(nc, f2_varid, 'units', 's-1');
netcdf.putAtt(nc, f2_varid, 'axis', 'Y');
netcdf.putAtt(nc, f2_varid, 'associates', 'frequency');


dir_varid=netcdf.defVar(nc, 'direction' ,'NC_FLOAT',[DIR]);
netcdf.putAtt(nc, dir_varid, 'long_name', 'sea surface wave to direction');
netcdf.putAtt(nc, dir_varid, 'standard_name', 'sea surface wave to direction');
netcdf.putAtt(nc, dir_varid, 'globwave_name', 'direction');
netcdf.putAtt(nc, dir_varid, 'units', 'degree');
netcdf.putAtt(nc, time_varid, 'axis', 'X');

EFTH_varid=netcdf.defVar(nc, 'efth' ,'NC_FLOAT',[DIR,FREQ,STAT,unlim_id]);
netcdf.putAtt(nc, EFTH_varid, 'units', 'm2 s rad-1');
netcdf.putAtt(nc, EFTH_varid, 'field', 'efth, scalar, series');
netcdf.putAtt(nc, EFTH_varid, 'long_name', 'sea surface wave directional variance spectral density');
netcdf.putAtt(nc, EFTH_varid, 'standard_name', 'sea_surface_wave_directional_variance_spectral_density');
netcdf.putAtt(nc, EFTH_varid, 'globwave_name', 'directional_variance_spectral_density');
netcdf.putAtt(nc, EFTH_varid, 'associates', 'time station frequency direction');
netcdf.putAtt(nc, EFTH_varid, 'content', 'TXYZ');


d_varid=netcdf.defVar(nc, 'dpt' ,'NC_FLOAT',[STAT,unlim_id]);
netcdf.putAtt(nc, d_varid, 'long_name', 'depth');
netcdf.putAtt(nc, d_varid, 'standard_name', 'depth');
netcdf.putAtt(nc, d_varid, 'globwave_name', 'depth');
netcdf.putAtt(nc, d_varid, 'content', 'TX');
netcdf.putAtt(nc, d_varid, 'associates', 'time station');
netcdf.putAtt(nc, d_varid, 'units', 'm');


wnd_varid=netcdf.defVar(nc, 'wnd' ,'NC_FLOAT',[STAT,unlim_id]);
netcdf.putAtt(nc, wnd_varid, 'units', 'm s-1');
netcdf.putAtt(nc, wnd_varid, 'long_name', 'wind speed at 10m');
netcdf.putAtt(nc, wnd_varid, 'standard_name', 'wind speed');
netcdf.putAtt(nc, wnd_varid, 'globwave_name', 'wind speed');
netcdf.putAtt(nc, wnd_varid, 'content', 'TX');
netcdf.putAtt(nc, wnd_varid, 'associates', 'time station');


wnddir_varid=netcdf.defVar(nc, 'wnddir' ,'NC_FLOAT',[STAT,unlim_id]);
netcdf.putAtt(nc, wnddir_varid, 'units', 'degree');
netcdf.putAtt(nc, wnddir_varid, 'long_name', 'wind direction');
netcdf.putAtt(nc, wnddir_varid, 'standard_name', 'wind_from_direction');
netcdf.putAtt(nc, wnddir_varid, 'globwave_name', 'wind_from_direction');
netcdf.putAtt(nc, wnddir_varid, 'content', 'TX');
netcdf.putAtt(nc, wnddir_varid, 'associates', 'time station');


cur_varid=netcdf.defVar(nc, 'cur' ,'NC_FLOAT',[STAT,unlim_id]);
netcdf.putAtt(nc, cur_varid, 'units', 'm s-1');
netcdf.putAtt(nc, cur_varid, 'long_name', 'sea water speed');
netcdf.putAtt(nc, cur_varid, 'standard_name', 'sea_water_speed');
netcdf.putAtt(nc, cur_varid, 'globwave_name', 'sea_water_speed');
netcdf.putAtt(nc, cur_varid, 'content', 'TX');
netcdf.putAtt(nc, cur_varid, 'associates', 'time station');


curdir_varid=netcdf.defVar(nc, 'curdir' ,'NC_FLOAT',[STAT,unlim_id]);
netcdf.putAtt(nc, curdir_varid, 'units', 'degree');
netcdf.putAtt(nc, curdir_varid, 'long_name', 'direction from of sea water velocity');
netcdf.putAtt(nc, curdir_varid, 'standard_name', 'direction_of_sea_water_velocity');
netcdf.putAtt(nc, curdir_varid, 'globwave_name', 'direction_of_sea_water_velocity');
netcdf.putAtt(nc, curdir_varid, 'content', 'TX');
netcdf.putAtt(nc, curdir_varid, 'associates', 'time station');
                     


netcdf.endDef(nc);

netcdf.putVar(nc, time_varid,0,2, time-datenum('01011990 000000','ddmmyyyy HHMMSS'));
netcdf.putVar(nc, station_varid, 84);
netcdf.putVar(nc, station16_varid, nan*(ones(16,1)));
netcdf.putVar(nc, station_name_varid,pointID);
netcdf.putVar(nc, lon_varid, Lon);
netcdf.putVar(nc, lat_varid, Lat);
netcdf.putVar(nc, f_varid, freq);
netcdf.putVar(nc, f1_varid, frequency1);
netcdf.putVar(nc, f2_varid, frequency2);
netcdf.putVar(nc, dir_varid, dir);
netcdf.putVar(nc, EFTH_varid, EFTH);
netcdf.putVar(nc, d_varid, dpt);
netcdf.putVar(nc, wnd_varid, wndspd);
netcdf.putVar(nc, wnddir_varid, wnddir);
netcdf.putVar(nc, cur_varid, curspd);
netcdf.putVar(nc, curdir_varid, curdir);




netcdf.close(nc);


%%

fileattrib(ncfile,'+w');
ncwriteatt(ncfile,'efth','scale_factor', 1);
ncwriteatt(ncfile,'efth','add_offset', 0);
ncwriteatt(ncfile,'efth','valid_min', 0);
ncwriteatt(ncfile,'efth','valid_max', 1e+20);
ncwriteatt(ncfile,'efth','_FillValue', int32(9.97e+36));

ncwriteatt(ncfile,'dpt','scale_factor', 1);
ncwriteatt(ncfile,'dpt','add_offset', 0);
ncwriteatt(ncfile,'dpt','valid_min', -100);
ncwriteatt(ncfile,'dpt','valid_max', 1e+04);
ncwriteatt(ncfile,'dpt','_FillValue', int32(9.97e+36));

ncwriteatt(ncfile,'y','scale_factor', 1);
ncwriteatt(ncfile,'y','add_offset', 0);
ncwriteatt(ncfile,'y','valid_min', -90);
ncwriteatt(ncfile,'y','valid_max', 180);
ncwriteatt(ncfile,'y','_FillValue', int32(9.97e+36));

ncwriteatt(ncfile,'x','scale_factor', 1);
ncwriteatt(ncfile,'x','add_offset', 0);
ncwriteatt(ncfile,'x','valid_min', -180);
ncwriteatt(ncfile,'x','valid_max', 360);
ncwriteatt(ncfile,'x','_FillValue', int32(9.97e+36));

ncwriteatt(ncfile,'direction','scale_factor', 1);
ncwriteatt(ncfile,'direction','add_offset', 0);
ncwriteatt(ncfile,'direction','valid_min', 0);
ncwriteatt(ncfile,'direction','valid_max', 360);
ncwriteatt(ncfile,'direction','_FillValue', int32(9.97e+36));

                          
ncwriteatt(ncfile,'frequency','scale_factor', 1);
ncwriteatt(ncfile,'frequency','add_offset', 0);
ncwriteatt(ncfile,'frequency','valid_min', 0);
ncwriteatt(ncfile,'frequency','valid_max', 10);
ncwriteatt(ncfile,'frequency','_FillValue', int32(9.97e+36));

    
ncwriteatt(ncfile,'frequency1','scale_factor', 1);
ncwriteatt(ncfile,'frequency1','add_offset', 0);
ncwriteatt(ncfile,'frequency1','valid_min', 0);
ncwriteatt(ncfile,'frequency1','valid_max', 10);
ncwriteatt(ncfile,'frequency1','_FillValue', int32(9.97e+36));


ncwriteatt(ncfile,'frequency2','scale_factor', 1);
ncwriteatt(ncfile,'frequency2','add_offset', 0);
ncwriteatt(ncfile,'frequency2','valid_min', 0);
ncwriteatt(ncfile,'frequency2','valid_max', 10);
ncwriteatt(ncfile,'frequency2','_FillValue', int32(9.97e+36));


ncwriteatt(ncfile,'wnd','scale_factor', 1);
ncwriteatt(ncfile,'wnd','add_offset', 0);
ncwriteatt(ncfile,'wnd','valid_min', 0);
ncwriteatt(ncfile,'wnd','valid_max', 100);
ncwriteatt(ncfile,'wnd','_FillValue', int32(9.97e+36));

ncwriteatt(ncfile,'wnddir','scale_factor', 1);
ncwriteatt(ncfile,'wnddir','add_offset', 0);
ncwriteatt(ncfile,'wnddir','valid_min', 0);
ncwriteatt(ncfile,'wnddir','valid_max', 360);
ncwriteatt(ncfile,'wnddir','_FillValue', int32(9.97e+36));


ncwriteatt(ncfile,'cur','scale_factor', 1);
ncwriteatt(ncfile,'cur','add_offset', 0);
ncwriteatt(ncfile,'cur','valid_min', 0);
ncwriteatt(ncfile,'cur','valid_max', 100);
ncwriteatt(ncfile,'cur','_FillValue', int32(9.97e+36));

ncwriteatt(ncfile,'curdir','scale_factor', 1);
ncwriteatt(ncfile,'curdir','add_offset', 0);
ncwriteatt(ncfile,'curdir','valid_min', 0);
ncwriteatt(ncfile,'curdir','valid_max', 360);
ncwriteatt(ncfile,'curdir','_FillValue', int32(9.97e+36));
ncwriteatt(ncfile,'station','_FillValue', int32(-2147483647));
ncwriteatt(ncfile,'string16','_FillValue', int32(-2147483647));
end