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
%----------------------------------------------------------%
SWDEN.buoy_name=cellstr(flipud(rot90(ncread(ncf,'station_name'))));
SWDEN.time=convert_time(ncf,'time');
SWDEN.latitude=ncread(ncf,'latitude');
SWDEN.longitude=ncread(ncf,'longitude');
SWDEN.f=double(ncread(ncf,'frequency'));
SWDEN.Dir=ncread(ncf,'direction');
SWDEN.DENS=ncread(ncf,'efth');
DeltaDir=abs(SWDEN.Dir(2)-SWDEN.Dir(1));
DENS=SWDEN.DENS;
DENS=[DENS; DENS(1,:,:,:)];
for k=1:length(SWDEN.buoy_name)
for ii=1:length(SWDEN.f)
    SPEC(ii,k,:) = abs(trapz(0:pi*DeltaDir/180:2*pi,DENS(:,ii,k,:)));
end
end
SWDEN.SPEC=SPEC;
for k=1:length(SWDEN.buoy_name)
 Hs(k,:)=4.004*sqrt(abs(trapz(2*pi*SWDEN.f,SPEC(:,k,:))))/sqrt(2*pi);
end
 SWDEN.Hs=Hs;

for k=1:length(SWDEN.buoy_name)
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


