function topo=BoxSmoothTopo(fl,k)
   
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