function SmoothCoastalBoundaries(shpfl, LambdaKilometers, flout, MinArea,MinPerimeter)
%
%function SmoothCoastalBoundaries(shpfl, LambdaKilometers, flout, MinArea,MinPerimeter)
%
% Script to take coastline shapefile and smooth and resample the curves at
% uniform spacing.  Longitude and latitude of coastline points are treated as 
% functions of  distance along the coast in the smoothing and resampling operations. 
%
% Islands that are both smaller than a threshold area (MinArea- sq km) and whos 
% perimeter is shorter than a threshold length (MinPerimeter - km) are excluded. 
% from the smoothed shapefile. Lambda, MinArea, and MinPerimeter are in km
% 
%
% Inputs:
%       shpfl: filename of shapefile containing coastline to be smoothed
%       LambdaKilometers: Length in kilometers to smooth and resample the
%       coastline to
%       flout: File to write the smoothed subsampled coastline data to
%       MinArea(optional): Area in square kilometers to filter islands (see below)
%       MinPerimeter(optional): Length in kilometers to filter islands (see below)
%
%       Regarding MinArea and MinPerimeter: Islands that are both smaller 
%       than a threshold area (MinArea- sq km) and whos perimeter is 
%       shorter than a threshold length (MinPerimeter - km) are excluded
%       from the output shapefile. This is intended to filter out islands that
%	are too small to resolve in the mesh being built.
%
%


LambdaMeters=LambdaKilometers*1000;

if nargin < 4,
    MinArea=(LambdaKilometers*2)^2;
end
if nargin < 5,
    MinPerimeter= 4*LambdaKilometers;
end

SetPath
isplot=0

deg2km=111.132954
deg2rad=pi/180

earth=referenceSphere('Earth')

S = shaperead(shpfl);
N=length(S);
isisland=zeros(N,1);
for k=1:N
    ns(k)=length(S(k).X(1:end-1));
    if and(  S(k).X(end-1)==S(k).X(1) , S(k).Y(end-1)==S(k).Y(1) )
        isisland(k)=1;
    end
end
if sum(isisland)==N
    disp(['All features in ',GlobalCoastlineFile,' are closed islands'])
else
    disp([int2str(sum(isisland)),'  features in ',GlobalCoastlineFile,...
        ' are closed islands. out of:',int2str(N),' total features'])
end

for k=1:N
    x=S(k).BoundingBox(:,1);
    y=S(k).BoundingBox(:,2);
    x0(k)=mean(x);
    y0(k)=mean(y);
    dx(k)=abs(x(1)-x(2));
    dy(k)=abs(y(1)-y(2));
end
A=[dx.*cos(deg2rad*y0/180)*deg2km].*[dy*deg2km];
k=find(A>.500*.500);
S=S(k);
ns=ns(k);
N=length(S);
[ss,k]=sort(-ns);
S=S(k);
ns=ns(k);
N=length(S);
n=0;

for k=1:N
    x=S(k).X(1:end-2);%unique* points
    y=S(k).Y(1:end-2);
%deal with "messy" OpenStreetMap coasts
    jb=find(isnan(x+y));
    if ~isempty(jb)
        x=x(1:[min(jb)-1]);
        y=y(1:[min(jb)-1]);%truncate coastline at first discontinuity
    end
    zz=x+i*y;
    [zzu,j]=unique(zz);%find unique points before end
    j=sort(j);
    if length(zzu)<length(zz),
        x=x(j);
        y=y(j);
    end

    if length(x)>2
        x=[x,x(1)];%close loop
        y=[y,y(1)];
        [xs,ys]=SmoothSubSampleCoastlineFast(x,y,LambdaMeters/10.,10);
        zs=xs+i*ys;
        dz=abs(zs(2:end)-zs(1:end-1));
        j=find(dz>1);
        if length(j)>0
            ji=1:j(1);
            xs=xs(ji);ys=ys(ji);
        end
    
        if length(xs)>2,
            n=n+1;
            S0(n).X=[xs(:);NaN]';
            S0(n).Y=[ys(:);NaN]';
            S0(n).X0=S(n).X;
            S0(n).Y0=S(n).Y;
            S0(n).area = areaint(ys,xs,earth) / 10^6;
            S0(n).perim=sum( deg2km*abs(  cos(deg2rad*mean(ys))*(  xs(2:end)-xs(1:end-1) ) + i*(ys(2:end)-ys(1:end-1)) ) );
            S0(n).Geometry=S(k).Geometry;
        end
    end

    if mod(k,1000)==0,
        disp(['progress a:', num2str(  sum(ns(1:k)) / sum(ns) ), ', progress b:', num2str(  k/N  )]);
    end
end

N=length(S0)
for k=1:N
    ns0(k)=length(S0(k).X(1:end-1));
end
S=S0;

eval(['save -v7.3 ',flout,'.mat']);

%filter out small islands
js=[ ];
jsp=[ ];
clear A P
for n=1:length(S),
    A(n)=S(n).area;
    if S(n).area>MinArea,
        js=[js,n];
    end
    %if isempty(S(n).perim),S(n).perim=0;end
    P(n)=S(n).perim;
    if S(n).perim>MinPerimeter,
        jsp=[jsp,n]; % perimeter criteria for narrow Atolls
    end
end

if isplot
    close all
    figure;
    for k=1:length(js)
        n=js(k);
        xs=S(n).X(1:end-1);
        ys=S(n).Y(1:end-1);
        plot(xs,ys,'k-');hold on
    end
    pindx=setdiff(jsp,js) % narrow attols
    for k=1:length(pindx)
        n=pindx(k);
        xs=S(n).X(1:end-1);
        ys=S(n).Y(1:end-1);
        plot(xs,ys,'r.-');
    end
end

jsg=union(js,jsp);

S=S(jsg);
eval(['save -v7.3 ',flout,'.mat']);
S=rmfield(S,'X0')
S=rmfield(S,'Y0')

shapewrite(S, [flout,'.shp']);

%Output in jigsaw .msh format
%BoundaryShape2msh(S,[flout,'.msh']);
