function [xi,yi]=interpCoastline(x,y,dsmooth,itz);
%function [xs,ys]=SmoothSubSampleCoastline(x,y,dsmooth,lambda);
%interpolate coastline segment at uniform dsmooth (m) distance. 
%   input:
%       x,y :   longitude, latitude points along a coastline. 
%               If the coastline defines an island 
%               then x(1)==x(end) and y(1)==y(end),
%       dsmooth: distance in meters to interpolate the coastline to
%       itz: replace points at longitude 180, -180 with nans to supress artifacts
%
%   outputs : 
%       xi,yi : longitude, latitude points the coastline resampled at dsmooth distance
%
if nargin<4,itz=0;end
if itz,j=find(or(x==-180,x==180));end
z=x+i*y;
if z(1)==z(end),
    isisland=1;
else
    isisland=0;
end
lat2m=110574.;
dx=x(2:end)-x(1:end-1);
dy=y(2:end)-y(1:end-1);
ymp=(y(2:end)+y(1:end-1))/2;
lon2m=111320.*cos(ymp*pi/180);
d=sqrt(  (dx.*lon2m).^2 + (dy.*lat2m).^2   ); 
d=[0,cumsum(d)];%distance
di=d(1):dsmooth:d(end);
if itz,z(j)=NaN;end
%zi=interp1(d,z,di);
[du,ju]=unique(d);
ju=sort(ju);
zi=interp1(d(ju),z(ju),di);%degenerate cases in OSM Antarctica
if isisland
    zi=[zi,zi(1)];
end
xi=real(zi);
yi=imag(zi);

