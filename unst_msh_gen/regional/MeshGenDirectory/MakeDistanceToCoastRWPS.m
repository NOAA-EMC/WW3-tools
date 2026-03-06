%This script creates a "Distance to Coast file" using a US coastline file 
% and custom specification of other designated points.  The PSLG intended for creation of the boundary
% is used as well as a global bathymetry file and US coastline shapefile.

SetPath
isplot=0

%Input files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PSLGfile='GlobalCoastlineGSHHS.PSLG.msh'
%PSLGfile='GlobalCoastlineOSM.PSLG.msh'
%PSLGfile='PSLGboundaryOSMxGSHHS.BOXES.msh'

%Coastline Boundary file
%PSLGfile='GlobalCoastlineGSHHS.PSLG.msh'
%Global Bathymetry file (netcdf)
%GlobalTopoFile='/scratch3/NCEPDEV/climate/Keston.Smith/RWPS/Data/RTopo_2_0_4_GEBCO_v2023_60sec_pixel.nc'
%Shape file with US coastline 
%TargetCoastlineFile='/scratch3/NCEPDEV/climate/Keston.Smith/RWPS/Data/us_coastline/tl_2023_us_coastline.shp'

%Output files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output "distance to coast" grided jigsaw .msh file
DFunOutFile=['DFun.',PSLGfile(1:end-4),'.msh']
%Output bathymetry on same grid as DFunOutFile
TopoFile=['Topo.',DFunOutFile];
%Temporary matlab file
FileOutMatlab=[PSLGfile(1:end-4),'.MakeDistance.mat']


US=shaperead(TargetCoastlineFile);
NUS=length(US)
xus=[];
yus=[];
x0=[]
y0=[];
interpDist=100. %meters
SmoothN=25;% box car smooth coastline at width 2(*SmoothN+1)*interpDist 
% then subsample soothed coastline every SmoothN points to create smaller
% coastline data set. i.e. remove small scale variation
for k=1:NUS
    x=US(k).X(1:end-1);
    y=US(k).Y(1:end-1);
    x0=[x0,x];
    y0=[y0,y];
%    xs=x;ys=y; %1.5 mil points
    if length(x)>3
        [xs,ys]=SmoothSubSampleCoastlineFast(x,y,interpDist,SmoothN);
        xus=[xus,xs];
        yus=[yus,ys];
    end
end
n0=length(xus);

[xpi,ypi]=ExtraPointsOfIntrest
 
xus=[xus(:)',xpi(:)']
yus=[yus(:)',ypi(:)']

Blon=[129.91 10.71];
Blat=[-30.42 79.99];
Blon=LonCon(Blon)


if isplot,
    xusp=LonCon(xus);
    clf;
    plot(xusp(1:n0)    ,yus(1:n0),'k.',...
     xusp(n0+1:n1) ,yus(n0+1:n1),'co', ...
     xusp(n1+1:end),yus(n1+1:end),'rx');
    hold on
    BoundingBox(Blon-360,Blat,'r');
    grid on;
    title('Updated RWPS high res target points');
    kprint('RWPSHighResPointsX.jpg');
end
topo=BoxSmoothTopo(GlobalTopoFile,2);

lon = topo.point.coord{:,1};
lat = topo.point.coord{:,2};

%p=loadmsh('PSLGboundaryOSMxGSHHS1km.msh');
p=loadmsh(PSLGfile);

%p.x=p.point.coord(:,1)-360;
p.x=p.point.coord(:,1);
p.y=p.point.coord(:,2);
p.x=LonCon(p.x);
p.edges=p.edge2.index(:,1:2);
plot(p.x,p.y,'g.');

%Remove near boundary points- Not land targets of resolution
dx=1;
j=find(p.x<max(p.x)-dx);p=subpslgFast(p,j);
j=find(p.x>min(p.x)+dx);p=subpslgFast(p,j);
j=find(p.y<max(p.y)-dx);p=subpslgFast(p,j);
j=find(p.y>min(p.y)+dx);p=subpslgFast(p,j);
plot(p.x,p.y,'y.');

%j=find(xus>90);
%xus(j)=xus(j)-360;

xus=LonCon(xus);
p.x=LonCon(p.x);
%This loop Takes a few min:Find PSLG points near us coastline
deg2km=111.132954
np=length(p.x);
NP=10000;

t0=now;
clear d;
nus=length(xus);
for k=1:NP:np
    j=k:min(np,k+NP);
    x=p.x(j);
    y=p.y(j);
    d(j)=min(deg2km*abs([x(:)-xus(:)'].*[cos(pi*y(:)/180)*ones(1,nus)]+i*[y(:)-yus(:)'])');
    %if mod(k,1000)==0,k/np,end
    if mod(k,1)==0,
        t1=now;
        est=(np-k)*[(t1-t0)/k];
        disp(['Computing PSLG distance to coast, etimated time remaing: ' ,num2str(est*24*60),' miniutes']);
    end
end

dmin=100;
%dmin=50;
j=find(d<dmin);
plot(p.x(j),p.y(j),'b.')

xUSsl=p.x(j);%"US" shore line points
yUSsl=p.y(j);

whos x xus xsl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Some ploting stuff

clf;
j=find(xus>90);
xus(j)=xus(j)-360;
Blon=[129.91-360, 10.71];
Blat=[-30.42, 79.99];
if isplot,
    plot(p.x,p.y,'k.');
    hold on
    BoundingBox(Blon,Blat,'k'); 
    th=title('PSLG defining boundary');
    set(th,'FontSize',22);
    kprint('PSLG.jpg')

    plot(xus,yus,'r.');
    th=title('PSLG with US coasline data points(red)');
    set(th,'FontSize',22);
    kprint('PSLGwUSpoints.jpg')

    j=find(d<100);
    plot(p.x(j),p.y(j),'c.');
    th=title('PSLG with US coasline data points and PSLG points near us (cyan)');
    set(th,'FontSize',22);
    kprint('PSLGwUSpointsPSLGnear.jpg')
end



xx=[xUSsl(:);xus(:)]';
yy=[yUSsl(:);yus(:)]';
xx=LonCon(xx);
plot(xx,yy,'c.');

save tmp.mat

D=DistanceToCoast(lon,lat,xx,yy);

eval(['save -v7.3 ',FileOutMatlab, ' lon lat D xx yy xus yus xUSsl yUSsl']) ;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make .msh files for jigsaw 

topo=BoxSmoothTopo(GlobalTopoFile,2);

Dfun=topo;
 
lon=Dfun.point.coord{:,1};
j=find(lon<90);lon(j)=lon(j)+360;
j0=setdiff(1:length(lon),j);j0=j0(:);
lon=lon([j0(:);j(:)]);
D1=[D(:,j0),D(:,j)];
Dfun.point.coord{:,1}=lon;
Dfun.value=D1;

Dmax=20004000%max distance between two points on earth
D1(find(D1>Dmax))=Dmax;
Dfun.value=D1;

figure;clf;pcolor(Dfun.point.coord{:,1},Dfun.point.coord{:,2},exp(-Dfun.value/320000));
shading interp;colorbar;colormap('jet')
savemsh(DFunOutFile,Dfun);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save smoothed topo at same resolution as Distance function 
topo=BoxSmoothTopo(GlobalTopoFile,2);

D=topo.value;
lon=topo.point.coord{:,1};
j=find(lon<90);lon(j)=lon(j)+360;
j0=setdiff(1:length(lon),j);j0=j0(:);
lon=lon([j0(:);j(:)]);
D1=[D(:,j0),D(:,j)];
topo.point.coord{:,1}=lon;
topo.value=D1;
 
figure;clf;pcolor(topo.point.coord{:,1},topo.point.coord{:,2},topo.value);
shading interp;colorbar;colormap('jet')
savemsh(TopoFile,topo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
