function  A = EleArea(lon,lat,e)
%crude compute element Area
X=lon(e);
Y=lat(e);
mY=mean(Y')';
lat2m=111111.;
lon2m=lat2m.*cos(pi*mY/180.);
D3=abs( lon2m.*(X(:,2)-X(:,1)) + lat2m*i*( Y(:,2)-Y(:,1) ));
D2=abs( lon2m.*(X(:,3)-X(:,1)) + lat2m*i*( Y(:,3)-Y(:,1) ));
D1=abs( lon2m.*(X(:,3)-X(:,2)) + lat2m*i*( Y(:,3)-Y(:,2) ));
S=(D1+D2+D3)/2;
A=sqrt(S .* (S - D1) .* (S - D2) .* (S - D3));
