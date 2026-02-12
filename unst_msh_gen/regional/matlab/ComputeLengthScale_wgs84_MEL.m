function lengthscale_wgs84_MEL=ComputeLengthScale_wgs84_MEL(x,y,e)
%function ComputeLengthScale(x,y,e)
% Input:
% x -lon [nn x 1]
% y -lat [nn x 1]
% e -element list[ ne x 3]
% 
% Output:
% lengthscale [ne x 1] lengthscale of element in km,
% defined by mean side length
mthd=1

x1=x(e(:,1));y1=y(e(:,1));
x2=x(e(:,2));y2=y(e(:,2));
x3=x(e(:,3));y3=y(e(:,3));

if mthd==1
    wgs84 = wgs84Ellipsoid("km");
    D3 = distance([y1,x1],[y2,x2],wgs84);
    D1 = distance([y2,x2],[y3,x3],wgs84);
    D2 = distance([y3,x3],[y1,x1],wgs84);
    lengthscale_wgs84_MEL=[D1+D2+D3]/3;% mean edge length
  %  lengthscale_wgs84_MEL=min(min(D1,D2),D3);% shortest edge
else
    R=6378.100 %radius of earth in km
    D03 = distance([y1,x1],[y2,x2],'degrees');
    D01 = distance([y2,x2],[y3,x3],'degrees');
    D02 = distance([y3,x3],[y1,x1],'degrees');
    D0=[D01+D02+D03]/3;
    lengthscale_wgs84_MEL=R*sin(D0*pi/180);
end