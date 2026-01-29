function MakeCoastalBoundariesGSHHS

% make coastlines for mesh generation.  
%Islands smaller than a threshold area are excluded.
%Narrow islands(some attols) with less than critical area are included if 
%perimeter is longer than minPerimeter.  Global and Pacific values are treated differently
%Coastline is smoothed and subsampled to usefull scale for mesh generation in RWPS
%various smoothings of coastlines


clear
isplot=0
deg2km=111.132954
deg2rad=pi/180


%Pacific 
minareaP=1;% square km
minPerimeterP=3;%  km

minareaG=1;% square km
minPerimeterG=10;% km

geom.mshID='EUCLIDEAN-MESH'
geom.fileV = 3
%filter out small islands
PacLon=[-140,140]
PacLat=[-20,40];

earth=referenceSphere('Earth')

gcfl='../RWPS/Data/GlobalCoast/GSHHS_shp/f/GSHHS_f_L1.shp'
S = shaperead(gcfl);
N=length(S);
isisland=zeros(N,1);
for k=1:N
    ns(k)=length(S(k).X(1:end-1));
    if and(  S(k).X(end-1)==S(k).X(1) , S(k).Y(end-1)==S(k).Y(1) )
        isisland(k)=1;
    end
end
if sum(isisland)==N
    disp(['All features in ',gcfl,' are closed islands'])
else
    disp([int2str(sum(isisland)),'  features in ',gcfl,...
        ' are closed islands. out of:',int2str(N),' total features'])
end

%S0=S;%???
n=0;
for k=1:N
    x=S(k).X(1:end-1);
    y=S(k).Y(1:end-1);
%    [xs,ys]=SmoothSubSampleCoastlineFast(x,y,100.,75);
    [xs,ys]=SmoothSubSampleCoastlineFast(x,y,50.,10);%500 m coastline
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
            
    if mod(k,1000)==0,
        disp(['progress a:', num2str(  sum(ns(1:k))/sum(ns)  )]);
        disp(['progress b:', num2str(  k/N  )]);
    end
end

N=length(S0)
for k=1:N
    ns0(k)=length(S0(k).X(1:end-1));
end
S=S0;

save -v7.3 GlobalCoastlineGSHHS.mat S ns0
%filter out small islands
js=[154,53];
jsp=[154,53];
clear A P
for n=1:length(S),
    xp=mean(S(n).X(1:end-1));
    yp=mean(S(n).Y(1:end-1));
    if and( or(xp<PacLon(1),xp>PacLon(2)), and(yp>PacLat(1),yp<PacLat(2)) )
        minarea=minareaP;
        minperim=minPerimeterP;
    else
        minarea=minareaG;
        minperim=minPerimeterG;
    end
    A(n)=S(n).area;
    if S(n).area>minarea,
        js=[js,n];
    end
    %if isempty(S(n).perim),S(n).perim=0;end
    P(n)=S(n).perim;
    if S(n).perim>minperim,
        jsp=[jsp,n]; % perimeter criteria for narrow Atolls
    end
    Amin(n)=minarea;
    Pmin(n)=minperim;
    if mod(n,1000)==0,
        disp(['progress ', num2str(  n/ length(S) )]);
     end

end

close all
if isplot
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

%js=find(A>minarea);
S=S(jsg);
save -v7.3 GlobalCoastlineGSHHS.mat S
S=rmfield(S,'X0')
S=rmfield(S,'Y0')

BoundaryShape2msh(S,'GlobalCoastlineGSHHS.msh');
shapewrite(S, 'GlobalCoastlineGSHHS.shp');
