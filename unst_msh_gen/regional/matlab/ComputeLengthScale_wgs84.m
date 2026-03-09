function lengthscale_wgs84=ComputeLengthScale_wgs84_area(x,y,e)
%function ComputeLengthScale(x,y,e)
% Input:
% x -lon [nn x 1]
% y -lat [nn x 1]
% e -element list[ ne x 3]
% 
% Output:
% lengthscale [ne x 1] lengthscale of element in m based on 
% element area assuming equilateral

EqTrC=4/sqrt(3);
[ne,three]=size(e)
lengthscale=zeros(1,ne);
wgs84 = wgs84Ellipsoid("km");
n=length(x);
xp=x;
yp=y;
xp(n+1)=NaN;
yp(n+1)=NaN;
[ne,three]=size(e)
eP=[e(:,1:3),n+1+zeros(ne,1)];
Xp=xp(eP)';
Yp=yp(eP)';
t0=now
clear a
a = areaint(Yp(:),Xp(:),wgs84);
t1=now;
60*24*(t1-t0)
et=60*24*(t1-t0);
disp([' total time : ',num2str(et),' min']);
lengthscale_wgs84=sqrt(a*EqTrC);
