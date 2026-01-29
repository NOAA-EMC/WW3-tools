function topo=topo2msh(fl,flout)

x=ncread(fl,'lon');
y=ncread(fl,'lat');
z=ncread(fl,'bed_elevation');
x=[x(2:end)+x(1:end-1)]/2;
y=[y(2:end)+y(1:end-1)]/2;

%z = reshape( z,length(y), length(x));
topo.point.coord{:,1}=x(:);
topo.point.coord{:,2}=y(:);
topo.value=z;

topo.mshID='ELLIPSOID-GRID'
topo.fileV=3;

nargin
if nargin>2
    savemsh(flout,topo);
end
 