function [SWDEN] = swden_ww3_read(ncf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function reads WW3 directional spectral density file %
% netcdf format                                             %
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
% f: frequency (Hz)
% Dir: direction (degree)
% DENS: Directional Spectral Density [freq,direction,stations,time]
% SPEC: Spectral Density [freq,stations,time]
% Hs: Significant Wave Heigth (m)
% Fp: Peak Freq (Hz)
% mwvdir: mean wave direction (degrees)
% cur: current speed (m/s)
% curdir: current direction (deg)
% wnd: wind speed (m/s)
% wnddir: wind direction (deg)
% dpt: depth (m)

%----------------------------------------------------------%
SWDEN.buoy_name=cellstr(flipud(rot90(ncread(ncf,'station_name'))));
SWDEN.time=convert_time(ncf,'time');
SWDEN.latitude=ncread(ncf,'latitude');
SWDEN.longitude=ncread(ncf,'longitude');
SWDEN.f=double(ncread(ncf,'frequency'));
SWDEN.Dir=ncread(ncf,'direction');
SWDEN.DENS=ncread(ncf,'efth');
SWDEN.cur=ncread(ncf,'cur');
SWDEN.curdir=ncread(ncf,'curdir');
SWDEN.wnd=ncread(ncf,'wnd');
SWDEN.wnddir=ncread(ncf,'wnddir');
SWDEN.dpt=ncread(ncf,'dpt');
DeltaDir=abs(SWDEN.Dir(2)-SWDEN.Dir(1));
DENS=SWDEN.DENS;
DENS=[DENS; DENS(1,:,:,:)];
[DD,FF,KK,TT]=size(SWDEN.DENS);

DENSCOS=cosd([SWDEN.Dir;SWDEN.Dir(1)]).*DENS;
DENSSIN=sind([SWDEN.Dir;SWDEN.Dir(1)]).*DENS;
for k=1:KK
for ii=1:FF
SPEC(ii,k,:) = abs(trapz(0:pi*DeltaDir/180:2*pi,DENS(:,ii,k,:)));
SPECA(ii,k,:) = (trapz(0:pi*DeltaDir/180:2*pi,DENSCOS(:,ii,k,:)));
SPECB(ii,k,:) = (trapz(0:pi*DeltaDir/180:2*pi,DENSSIN(:,ii,k,:)));
end
end
SWDEN.SPEC=SPEC;
for k=1:KK
 Hs(k,:)=4.004*sqrt(abs(trapz(2*pi*SWDEN.f,SPEC(:,k,:))))/sqrt(2*pi);
 AA(k,:)=((trapz(2*pi*SWDEN.f,SPECA(:,k,:))));
 BB(k,:)=((trapz(2*pi*SWDEN.f,SPECB(:,k,:))));
end
 SWDEN.Hs=Hs;
 SWDEN.mwvdir=(180*atan2(BB,AA)/pi)+180;
for k=1:KK
    clear SPECmax
    clear SPECtmp
    SPECtmp(:,:)=SPEC(:,k,:);
    SPECmax=nanmax(SPECtmp);
for i=1:length(SWDEN.time)
    clear i1
    clear j1
    [i1,j1]=find(SPECtmp(:,i)==SPECmax(i));
    Fp(k,i)=nanmean(SWDEN.f(i1));
end
SWDEN.Fp=Fp;    
end


