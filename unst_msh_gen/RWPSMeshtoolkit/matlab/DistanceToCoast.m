function D=DistanceToCoast(lon,lat,lonP,latP)
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
LambdaBox=10;%only use points within 5 degrees lat for 10x+ speed up
latP=latP(:);
lonP=lonP(:);

lonP=LonCon(lonP);
lon=LonCon(lon);

for j=1:ny % this loop takes ~16 hours
   t00=now;
   j/ny
    lon2m=single(111320.*cos(lat(j)*pi/180));
    jbox=find(abs(lat(j)-latP)<LambdaBox);
    if ~isempty(jbox)
        DLON= mod(lon(:)'-lonP(jbox),360) ;% large matrix
        %DLON= mod(lon(:)'-lonP(:),360) ;% large matrix
        %ld1=min (  abs(  [ DLON ]*lon2m + i*[ones(1,nx)*lat(j) - latP(:)]*lat2m  ) );
        %ld2=min (  abs(  [360 - DLON ]*lon2m + i*[ones(1,nx)*lat(j) - latP(:)]*lat2m  ) ); 360-ld1;
        ld1=min (  abs(  [ DLON ]*lon2m + i*[ones(1,nx)*lat(j) - latP(jbox)]*lat2m  ) );
        ld2=min (  abs(  [360 - DLON ]*lon2m + i*[ones(1,nx)*lat(j) - latP(jbox)]*lat2m  ) ); 
        D(j,:)=min(  [ ld1(:),ld2(:) ]'  ); % take shortes distance of east or west distance
    else
        D(j,:)=inf+zeros(1,nx);
    end
    t1=now;
    est=(ny-j)*[(t1-t00)];
    if mod(j,1)==0,
        disp(['Processing for latitude: ',num2str(lat(j)),' ,  data points: ',int2str(length(jbox)) ]);
        disp(['etimated 0 time remaing: ' ,num2str(est*24),' hrs']);
        est=(ny-j)*[(t1-t0)/j];
        disp(['etimated time remaing: ' ,num2str(est*24),' hrs']);
    end
end

