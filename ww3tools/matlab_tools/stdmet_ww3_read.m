function [STDMET] = stdmet_ww3_read(ncf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function reads WW3 wave bulk stats time series       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Ali Abdolali Feb 2023 ali.abdolali@noaa.gov          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  input data %--------------------------------------------%
% ncf: name of netcdf file
%  output data %--------------------------------------------%
% buoy_name: list of buoy
% time (Matlab time)
% longitude of points (Degree)
% latitude of points (Degree)
% mwvdir: mean wave direction (deg)
% Hs: Significant Wave Heigth (m)
% Fp: Peak Freq (Hz)
%----------------------------------------------------------%
STDMET.buoy_name=cellstr(flipud(rot90(ncread(ncf,'station_name'))));
STDMET.time=convert_time(ncf,'time');
STDMET.latitude=ncread(ncf,'latitude');
STDMET.longitude=ncread(ncf,'longitude');
STDMET.Hs=ncread(ncf,'hs');
STDMET.Fp=ncread(ncf,'fp');
STDMET.mwvdir=ncread(ncf,'th1m');




