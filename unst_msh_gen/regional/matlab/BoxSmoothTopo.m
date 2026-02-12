function topo=BoxSmoothTopo(fl,k)
%function topo=BoxSmoothTopo(fl,k)
% Apply 2D boxcar smoothing and sub sampling of bathymetry in file fl.
% jigsaw can sometimes run into memory problems with large stuctures 
% of this type when generating non-global meshes. Retuns topo, a 
% gridded jigsaw structure.
%   inputs: 
%       fl : filename pointing to a netcdf bathymetry file with variables
%           lon : length (nx+1) longitude of grid
%           lat : length (ny+1) latitude of grid
%           bed_elevation : ( nx by ny)  average bathymetric depth in cell
%       k : integer to resample the bathymetry at.The smoothing is carried
%           out using a 2 D Boxcar smoother (square average) with sides of
%           length (2 k + 1) and resampled every k points
%   output:
%       topo :  jigsaw format girdded data representing the smoothed resampled 
%               data.  
%

 
lon=ncread(fl,'lon');
lat=ncread(fl,'lat');
z=ncread(fl,'bed_elevation');

lon=[lon(1:end-1)+lon(2:end)]/2;
lat=[lat(1:end-1)+lat(2:end)]/2;
if k==0,
    topo.point.coord{:,1}=lon;
    topo.point.coord{:,2}=lat;%fix to centers later
    topo.value=double(z');
    
    topo.mshID='ELLIPSOID-GRID'
    topo.fileV=3;

else
    
    n=2*k+1;       
    w=ones(n,n)/n/n;
    [nx,ny] = size(z)
    zw=conv2(z,w,'same');
    zw=zw(k+1:n:nx-k,k+1:n:ny-k);
    lonw=lon(k+1:n:nx-k);
    latw=lat(k+1:n:ny-k);
    
    topo.point.coord{:,1}=lonw;
    topo.point.coord{:,2}=latw;%fix to centers later
    topo.value=double(zw');
    
    topo.mshID='ELLIPSOID-GRID'
    topo.fileV=3;
end
