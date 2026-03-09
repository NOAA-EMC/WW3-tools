function topo=topo2msh(fl,flout)
% convert bathymetry in netcdf file fl to jigsaw m.msh format
% Retuns topo, a gridded jigsaw structure.
%   inputs: 
%       fl : filename pointing to a netcdf bathymetry file with variables
%           lon : length (nx+1) longitude of grid
%           lat : length (ny+1) latitude of grid
%           bed_elevation : ( nx by ny)  average bathymetric depth in cell
%       flout : jigsaw .msh file to write topo to.
%   output:
%       topo :  jigsaw format girdded data representing the data in fl
%

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
 
