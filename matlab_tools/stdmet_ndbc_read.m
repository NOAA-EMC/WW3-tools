function [STDMET] = stdmet_ndbc_read(ncf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function reads NDBC standard Met time series         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Ali Abdolali Feb 2023 ali.abdolali@noaa.gov          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  input data %--------------------------------------------%
% ncf: name of netcdf file
%  output data %--------------------------------------------%
% time (Matlab time)
% longitude of points (Degree)
% latitude of points (Degree)
% windspd: wind speed
% winddir: wind direction
% gust: wind gust
% mwvdir: mean wave direction (deg)
% Hs: Significant Wave Heigth (m)
% Tp: Peak Period (s)
% Tm: Mean Period (s)
% P: pressure at MSL
% AT: air temperature
% SST: sea surface temperature
% dewptT: dewpoint temperature
% visibility: visibility
% elev: water surface elevation
%----------------------------------------------------------%
STDMET.time=convert_time(ncf,'time');
STDMET.latitude=ncread(ncf,'latitude');
STDMET.longitude=ncread(ncf,'longitude');
windspd=ncread(ncf,'wind_spd');
winddir=ncread(ncf,'wind_dir');
gust=ncread(ncf,'gust');
Hs=ncread(ncf,'wave_height');
Tp=ncread(ncf,'dominant_wpd');
Tm=ncread(ncf,'average_wpd');
mwvdir=ncread(ncf,'mean_wave_dir');
P=ncread(ncf,'air_pressure');
AT=ncread(ncf,'air_temperature');
SST=ncread(ncf,'sea_surface_temperature');
dewptT=ncread(ncf,'dewpt_temperature'); 
visibility=ncread(ncf,'visibility');
elev=ncread(ncf,'water_level'); 
%----------------------------------------------------------%
STDMET.windspd(:,1)=windspd(1,1,:);
STDMET.winddir(:,1)=winddir(1,1,:);
STDMET.gust(:,1)=gust(1,1,:);
STDMET.Hs(:,1)=Hs(1,1,:);
STDMET.Tp(:,1)=Tp(1,1,:);
STDMET.Tm(:,1)=Tm(1,1,:);
STDMET.mwvdir(:,1)=mwvdir(1,1,:);
STDMET.P(:,1)=P(1,1,:);
STDMET.AT(:,1)=AT(1,1,:);
STDMET.SST(:,1)=SST(1,1,:);
STDMET.dewptT(:,1)=dewptT(1,1,:); 
STDMET.visibility(:,1)=visibility(1,1,:);
STDMET.elev(:,1)=elev(1,1,:); 


