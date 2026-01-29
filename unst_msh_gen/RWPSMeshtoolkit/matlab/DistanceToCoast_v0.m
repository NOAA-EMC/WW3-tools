function D=DistanceToCoast_v0(lon,lat,lonP,latP)
%function D=DistanceToCoast(lon,lat,lonP,latP) 
% Computes grided distance to pointset on Globe
%
% Computes shortest distance to pointset [lonP,latP] from coordinates of 
% lat by lon grid in m. 
%
% inputs:
%   lon -   [nx,1] lon coordinates for output matrix
%   lat -   [ny,1] lon coordinates for output matrix
%   lonP  - [nn,1] longitude coordiantes of point set
%   latP  - [nn,1] latitude coordiantes of point set
%
% outputs:
%   D   - [ny, nx] units meters
%   D(j,k)= minimumt Distance from point lat(j), lon(k) to pointset (lonP,latP)
%
clear d0
lat2m=single(110574.)
lon=single(lon);
lat=single(lat);
nx=length(lon)
ny=length(lat)
clear ld1 ld2 D;
t0=now;
latP=latP(:);
lonP=lonP(:);
for j=1:ny % this loop takes ~16 hours
    j/ny
    lon2m=single(111320.*cos(lat(j)*pi/180));
    DLON= mod(lon(:)'-lonP(:),360) ;% large matrix
    ld1=min (  abs(  [ DLON ]*lon2m + i*[ones(1,nx)*lat(j) - latP(:)]*lat2m  ) );
    ld2=min (  abs(  [360 - DLON ]*lon2m + i*[ones(1,nx)*lat(j) - latP(:)]*lat2m  ) ); 360-ld1;
    D(j,:)=min(  [ ld1(:),ld2(:) ]'  ); % take shortes distance of east or west distance
    t1=now;
    est=(ny-j)*[(t1-t0)/j];
    disp(['etimated time remaing: ' ,num2str(est*24),' hrs']);
end

