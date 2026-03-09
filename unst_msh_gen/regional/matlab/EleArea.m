function  A = EleArea(lon,lat,e)
%function  A = EleArea(lon,lat,e)
%Approximate element Area
%
%   inputs:
%       lon : (nn x 1) longitude coordinates of mesh nodes
%       lat : (nn x 1) latitude coordinates of mesh nodes
%       e   : (ne x 3) element matrix
%
%   output:
%       A : (ne x 1) aproximate area of each element in m^2
%
X=lon(e);
Y=lat(e);
mY=mean(Y')';
lat2m=110574.;
lon2m=111320.*cos(pi*mY/180.);
D3=abs( lon2m.*(X(:,2)-X(:,1)) + lat2m*i*( Y(:,2)-Y(:,1) ));
D2=abs( lon2m.*(X(:,3)-X(:,1)) + lat2m*i*( Y(:,3)-Y(:,1) ));
D1=abs( lon2m.*(X(:,3)-X(:,2)) + lat2m*i*( Y(:,3)-Y(:,2) ));
S=(D1+D2+D3)/2;
A=sqrt(S .* (S - D1) .* (S - D2) .* (S - D3));
