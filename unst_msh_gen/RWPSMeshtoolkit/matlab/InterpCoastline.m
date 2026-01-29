function [xi,yi]=interpCoastline(x,y,dsmooth,itz);
%function [xs,ys]=SmoothSubSampleCoastline(x,y,dsmooth,lambda);
%Sub sample coastline at dsmooth (m) distance. Boxcar smooth coastline to 
% dsmooth*lambda distance   
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
zi=interp1(d,z,di);
if isisland
    zi=[zi,zi(1)];
end
xi=real(zi);
yi=imag(zi);

