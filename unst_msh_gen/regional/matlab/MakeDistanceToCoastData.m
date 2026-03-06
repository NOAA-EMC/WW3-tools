%This script creates a "Distance to Coast file" using a US coastline file 
% and custom specification of other designated points.  The PSLG intended for creation of the boundary
% is used as well as a global bathymetry file and US coastline shapefile.

function MakeDistanceToCoastData(DX,TargetShape,xtrgt,ytrgt)

SetPath
if nargin<2
    TargetShape=TargetCoastlineFile
end
if nargin<4,
    xtrgt=[];
    ytrgt=[];
end

close all
isplot=0

%Input files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Coastline Boundary file
PSLGfile
GlobalTopoFile, %set as a global variable in SetPath
TargetCoastlineFile, %set as a global variable in SetPath

%Output files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output "distance to coast" grided jigsaw .msh file
DFunOutFile=['DFun.',PSLGfile(1:end-4),'.msh']
%Output bathymetry on same grid as DFunOutFile
TopoOutFile=['Topo.',DFunOutFile];
%Temporary matlab file
FileOutMatlab=[PSLGfile(1:end-4),'.MakeDistance.mat']

US=shaperead(TargetShape);
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

whos x0 xus

xus=[xus(:)',xtrgt(:)'];
yus=[yus(:)',ytrgt(:)'];

n0=length(xus);

p=loadmsh(PSLGfile);
p.x=p.point.coord(:,1);
p.x=LonCon(p.x);

p.y=p.point.coord(:,2);
p.edges=p.edge2.index(:,1:2);

xg=min(p.x-DX):DX:max(p.x+DX);
yg=min(p.y-DX):DX:max(p.y+DX);
nx=length(xg);
ny=length(yg);
XG=ones(ny,1)*[xg(:)'];
YG=yg(:)*ones(1,nx);

topo=BoxSmoothTopo(GlobalTopoFile,0);
lon = topo.point.coord{:,1};
lat = topo.point.coord{:,2};

localtopo=topo;
localtopo.point.coord{:,1}=xg(:);
localtopo.point.coord{:,2}=yg(:);
localtopo.value=interp2(topo.point.coord{:,1},topo.point.coord{:,2},topo.value,XG,YG);

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

%dmin=10;
dmin=50;
j=find(d<dmin);
plot(p.x(j),p.y(j),'y.')

xp=p.x(j);%"US" shore line points
yp=p.y(j);

xx=[xus(:)',xp(:)'];
yy=[yus(:)',yp(:)'];
D=DistanceToCoast(xg,yg,xx,yy);
 
save tmp.mat

savemsh(TopoOutFile,localtopo);

figure;clf;pcolor(localtopo.point.coord{:,1},localtopo.point.coord{:,2},localtopo.value);
shading interp;colorbar;colormap('jet');title('local topo')

localtopo.value=D;
savemsh(DFunOutFile,localtopo);

figure;clf;pcolor(localtopo.point.coord{:,1},localtopo.point.coord{:,2},exp(-localtopo.value/320000));
shading interp;colorbar;colormap('jet');title('Exp[- distance to target cost ]')
