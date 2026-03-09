%This script creates a "Distance to Coast file" using a US coastline file 
% and custom specification of other designated points.  The PSLG intended for creation of the boundary
% is used as well as a global bathymetry file and US coastline shapefile.

clear
close all
isplot=0

%Input files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PSLGfile='GlobalCoastlineGSHHS.PSLG.msh'
%PSLGfile='GlobalCoastlineOSM.PSLG.msh'
%PSLGfile='PSLGboundaryOSMxGSHHS.BOXES.msh'

%Coastline Boundary file
PSLGfile='GlobalCoastlineGSHHS.PSLG.msh'
%Global Bathymetry file (netcdf)
GlobalTopoFile='/scratch3/NCEPDEV/climate/Keston.Smith/RWPS/Data/RTopo_2_0_4_GEBCO_v2023_60sec_pixel.nc'
%Shape file with US coastline 
TargetCoastlineFile='/scratch3/NCEPDEV/climate/Keston.Smith/RWPS/Data/us_coastline/tl_2023_us_coastline.shp'

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
%Now add Pacific Territories and COFA points
%------------------------------------------------------------------------
%NWS Pacific Region, via the Compact of Free Association, 
% %oversees operations of 5 Weather Service Offices (WSO) across 
% Micronesia in the Republic of Palau, Federated States of Micronesia 
% and the Republic of the Marshall Islands. WFO Guam provides routine
% forecasts as well as WWA services for these areas. -Eric Lau
%Palau,Yap, Chuuk, Pohnpei, Majuro, Pago Pago, Wake island

%Wake Island
coord=[19.2796,166.6499];%E
xus=[xus,coord(2)];yus=[yus,coord(1)];
%Palau
coord=[7.4942, 134.5690]
xus=[xus,coord(2)];yus=[yus,coord(1)];
%Yap:9.5557° N, 138.1399° E
coord=[9.5557, 138.1399]
xus=[xus,coord(2)];yus=[yus,coord(1)];
%Chuuk: 7°25′N 151°47′E
coord=[7.374227, 151.754606]%E
xus=[xus,coord(2)];yus=[yus,coord(1)];
%Pohnpei › Coordinates: 6.8519° N, 158.2147° E
coord=[6.8519, 158.2147]
xus=[xus,coord(2)];yus=[yus,coord(1)];
%Majuro › Coordinates7.0667° N, 171.2667° E
coord=[7.0667, 171.2667]
xus=[xus,coord(2)];yus=[yus,coord(1)];
%Pago Pago › Coordinates: 14.2732° S, 170.7030° W
coord=[-14.2732, -170.7030]% SW
xus=[xus,coord(2)];yus=[yus,coord(1)];
% Kosrae › Coordinates :5.3096° N, 162.9815° E
coord=[5.3096, 162.9815]
xus=[xus,coord(2)];yus=[yus,coord(1)];
%------------------------------------------------------------------------

%Marshall Islands
%Majuro
coord=[7.0667, 171.2667]
%xus=[xus,coord(2)];yus=[yus,coord(1)];
%Ebeye
coord=[8.7815, 167.7373]
%xus=[xus,coord(2)];yus=[yus,coord(1)];
%Micronesia
%Kolonia,
coord=[6.9636, 158.2102]
%xus=[xus,coord(2)];yus=[yus,coord(1)];
%Pohnpei,
coord=[6.8519, 158.2147]
%xus=[xus,coord(2)];yus=[yus,coord(1)];
% Chuuk-Weno
coord=[7.4523, 151.8422]
%xus=[xus,coord(2)];yus=[yus,coord(1)];
% Tofol
coord=[5.3256, 163.0086]
%xus=[xus,coord(2)];yus=[yus,coord(1)];
%Colonia -between Palau and Guam
coord=[9.5164,138.1222]
%xus=[xus,coord(2)];yus=[yus,coord(1)];
n1=length(xus);
%------------------------------------------------------------------------
%Points from Curt
%Howland Island -Baker
coord=[0.8113, -176.6183]%W
xus=[xus,coord(2)];yus=[yus,coord(1)];
%Johnston atoll
coord=[16.7295, -169.5336]%W
xus=[xus,coord(2)];yus=[yus,coord(1)];
%Palmyra
coord=[5.8885, -162.0787]%W
xus=[xus,coord(2)];yus=[yus,coord(1)];
%Jarvis Island
coord=[0.3744, -159.9967]%W
xus=[xus,coord(2)];yus=[yus,coord(1)];
%Baker Island
%0.1936° N, 176.4769° W
coord=[0.1936,-176.4769 ]%W
xus=[xus,coord(2)];yus=[yus,coord(1)];
%------------------------------------------------------------------------
 coord=[18.4101, -75.0115 ]%W
 xus=[xus,coord(2)];yus=[yus,coord(1)];
 %If your team is going to the effort to add resolution for Majuro in RMI, I would strongly suggest doing the same for these atolls in RMI:
 %1. Kwajalein Atoll, which is home to part of the Ronald Regard Ballistic Missile Test Site (https://home.army.mil/kwajalein/index.php) and the underprediction of wave heights by NWS' WaveWatchIII model in January 2024, which I discussed on our last call, is what led to significant damage at the base, per: https://www.youtube.com/shorts/jH-pGoQDdcg   
 %2. Enewetak Atoll, which is the home of the Runit Dome (https://en.wikipedia.org/wiki/Runit_Island), is threatened by wave-driven overwash that has serious implications for the US Department of State, Department of Defense/War, and Department of the Interior via the Intergovernmental Compact of Free Association (https://www.doi.gov/oia/compacts-of-free-association).
 %3. Bikini Atoll, for similar reasons as Enewetak, although the radionuclides are all over the place and not all dumped in one location.
 %I would note that most of the atolls drop off at a 70-80 degree slope from approximately 30 m depth (which is generally less than 1000 m from shore) to over 1000 m depth, so there is not a need for a large region of increasing resolution to capture a broad continental shelf as characterizes CONUS.
 
 %Kwajalein Atoll   - 9.1898° N, 167.4243° E
 %Enewetak Atoll     - 11.4654° N, 162.1890° E
 %Bikini Atoll       - 11.6065° N, 165.3768° E
 
 %Kwajalein Atoll   - 9.1898° N, 167.4243° E
 coord=[9.1898, 167.4243 ]%W
 xus=[xus,coord(2)];yus=[yus,coord(1)];
 %Enewetak Atoll     - 11.4654° N, 162.1890° E
 coord=[ 11.4654, 162.1890]%W
 xus=[xus,coord(2)];yus=[yus,coord(1)];
 %Bikini Atoll       - 11.6065° N, 165.3768° E
 coord=[11.6065, 165.3768 ]%W
 xus=[xus,coord(2)];yus=[yus,coord(1)];
 
Blon=[129.91 10.71];
Blat=[-30.42 79.99];
lon=Blon;j=find(lon<90);lon(j)=lon(j)+360;
Blon=lon;
xusp=xus;j=find(xus>0);xusp(j)=xus(j)-360;

if isplot,
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

p.x=p.point.coord(:,1)-360;
p.y=p.point.coord(:,2);
p.edges=p.edge2.index(:,1:2);
plot(p.x,p.y,'g.');

%Remove near boundary points- Not land targets of resolution
dx=1;
j=find(p.x<max(p.x)-dx);p=subpslgFast(p,j);
j=find(p.x>min(p.x)+dx);p=subpslgFast(p,j);
j=find(p.y<max(p.y)-dx);p=subpslgFast(p,j);
j=find(p.y>min(p.y)+dx);p=subpslgFast(p,j);
plot(p.x,p.y,'y.');

j=find(xus>90);
xus(j)=xus(j)-360;

%This loop Takes a few min:Find PSLG points near us coastline
deg2km=111.132954
np=length(p.x);
NP=10000;

clear d;
nus=length(xus);
for k=1:NP:np
    j=k:min(np,k+NP);
    x=p.x(j);
    y=p.y(j);
    d(j)=min(deg2km*abs([x(:)-xus(:)'].*[cos(pi*y(:)/180)*ones(1,nus)]+i*[y(:)-yus(:)'])');
    %if mod(k,1000)==0,k/np,end
    k/np
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%add points to off shore banks we want to refine here we have Georges Bank
%and banks around the Bahamas.

xFB= [ -79.9909  -78.1963  -78.3957]
yFB=  [23.7978   26.8928   24.1530]
xGB = -67.4517
yGB =   41.3567

xx=[xUSsl(:);xus(:);xFB(:);xGB(:)]';
yy=[yUSsl(:);yus(:);yFB(:);yGB(:)]';
j=find(xx>90);
xx(j)=xx(j)-360;
plot(xx,yy,'c.');

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
