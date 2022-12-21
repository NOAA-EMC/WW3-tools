function [dates]=convert_time(ncf,varid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program converts the time variable in the netcdf file relative to  %
% unit origin (unit attribute)to matlab time                              %
% Ali Abdolali (EMC/NCEP/NOAA ali.abdolali@noaa.gov                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%    INPUT    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ncf: the name of netcdf file 
%varid: the name of time variable 
%%%%%%%%%%%%%%%%%%%    OUTPUT    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dates: matlab time
%%%%%%%%%%%%%%%%%%%    example   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[matlab_time]=convert_time('ww3.nc','time');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dates     = [];
 times     = double(ncread(ncf,varid));
 timeunits =  ncreadatt(ncf,varid,'units');
 if (~ isempty(timeunits))
    if strcmp(timeunits,'days')
        dates = times;
    elseif strncmp(timeunits,'days since',numel('days since'))
        timebase   = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
         timeorigin = double(datenum(timebase(1),timebase(2),timebase(3)));
        dates      = times + timeorigin;
    elseif strncmp(timeunits,'hours since',numel('hours since'))
        timebase   = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
         timeorigin = double(datenum(timebase(1),timebase(2),timebase(3)));
        dates      = double(times ./ 24.0 + timeorigin);

    elseif strncmp(timeunits,'seconds since',numel('seconds since'))
        timebase   = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
        timeorigin = double(datenum(timebase(1),timebase(2),timebase(3)));
        dates      = double(times ./ 86400.0 + timeorigin);
     end
 end
end
