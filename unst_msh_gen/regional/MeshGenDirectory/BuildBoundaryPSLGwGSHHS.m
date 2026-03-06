% Script to build PSLG based on Global Self-consistent, Hierarchical, High-resolution Geography Database (GSHHG).
% Data available at https://www.ngdc.noaa.gov/mgg/shorelines/shorelines.html
% Before running this script run script MakeCoastalBoundariesGSHHS.m to smooth curves in GSHHS data to reflect 
% desired coastal resolution in final mesh.

CoastLineFile = 'GlobalCoastlineGSHHS.shp'
lonWest=129.91;lonEast=10.71;latSouth=-30.42;latNorth=79.99;

S=shaperead(CoastLineFile)
FileOutJigsaw=[CoastLineFile(1:end-4),'.PSLG.msh']
FileOutMatlab=[CoastLineFile(1:end-4),'.PSLGtmp.mat']
isplot=0;

%load GlobalCoastlineOSM.mat


%The overall objective of this project is to develop and implement into operations
% a  Regional Wave Prediction System (RWPS) that fulfills the needs of NWS marine coastal, 
% offshore, and high seas areas of responsibility. The domain of the RWPS will cover the Atlantic,
% Pacific, and Arctic oceans equal to the bounds of the Oceanic Domain as developed for the 
% National Blend of Models (NBM). For reference the NBM Oceanic Domain has corner points of LL - 30.42S - 129.91E UR - 79.99N - 10.71E. 

%Blon=[129.91 10.71];
%Blat=[-30.42 79.99];
Blon=[lonWest, lonEast];
Blat=[latSouth,latNorth];


N=length(S);
N1=N;
for k=1:N
    if mod(k,1000)==0,k/N,end
    lon=S(k).X(1:end-1);%remove trailing nan 
    lat=S(k).Y(1:end-1);
    ns(k)=length(lon);
    jWest=find(lon<90);
    if(length(jWest)==ns(k));
        S(k).X=[lon+360,NaN];
    end
    if and(0<length(jWest),ns(k)>length(jWest)),%make translated copy
        display(['duplicating coastal segment: ',int2str(k), ', nseg=',int2str(ns(k))])
        N1=N1+1;
        S(N1).X=[lon+360,NaN];
        S(N1).Y=[lat,NaN];
        ns(N1)=length(lon);
    end
end
N=length(S);
[tmp,j]=sort(-ns);
S=S(j);% sort to descending in length
ns=ns(j);% sort to descending in length

lon=Blon;j=find(lon<90);lon(j)=lon(j)+360;
Blon=lon;

%Make Bounding rectangle
Bx=[Blon(1),Blon(2),Blon(2),Blon(1),Blon(1)]
By=[Blat(1),Blat(1),Blat(2),Blat(2),Blat(1)]
             
clear pslg
sxp=[];syp=[ ];
xc=[];yc=[];
N=length(S);
pslg.x=[];
pslg.y=[];
pslg.edges=[];
nc=0;


minedeges=4 % minimum number of edges in land bounadry feature to include
minarea=1 % minimum area for islands to be included in boundary (not used in this version)

%Loop below finds intersection of boanary rectangle with coastline segments to create outer boandary


earth=referenceSphere('Earth');
DXM=10^10;%DXM=1
mdx=inf+ones(N,1);
for k=1:N
    isinbox(k)=0;
    x=S(k).X(1:end-1);% remove trailing nan (-1) and endpoint==startpoint (-2)
    y=S(k).Y(1:end-1);
    mx=mean(x);
    my=mean(y);
    if length(x)>minedeges,
        ji=find( insidepoly( x,y,Bx,By  ) );
        if length(ji)>2  
            [xi, yi,ii] = polyxpoly(x, y, Bx, By);
            xc=[xc;xi(:)];
            yc=[yc;yi(:)];
            if ~isempty(xi) %modify ob
                sxp=[sxp,x(1)];syp=[syp,y(1)];
                [iia,jja]=sort(ii(:,1));%sort to ascending order along segments
                xia=xi(jja);
                yia=yi(jja);
                iis=ii(jja,1);%sorted into ascending order along segments
                if ~insidepoly( x(1),y(1),Bx,By )% segment origonates outside box 
                    nseg=length(xia);
                    for j=1:2:nseg-1,
                        xs=[xia(j),x( iis(j)+1:iis(j+1) ),xia(j+1)];
                        ys=[yia(j),y( iis(j)+1:iis(j+1) ),yia(j+1)];
                        dx=abs(xs(2:end)-xs(1:end-1) +i*[ys(2:end)-ys(1:end-1)] );
                        mdx(k)=max(dx);

                        if max(dx)<DXM,
                            n0=length(pslg.x);
                            pslg.x=[pslg.x,xs];
                            pslg.y=[pslg.y,ys];
                            n1=length(pslg.x);
                            nedges=[[n0+1:n1-1];[n0+2:n1]]';
                            pslg.edges=[pslg.edges;nedges];
                            nc=nc+1;pslg.chains(nc).nodes=[n0+1:n1];
                            pslg.chains(nc).index=k;
                            pslg.chains(nc).type='starts outside';
                            pslg.chains(nc).BI=1;
                        end
                    end
                else % segment starts inside outer boundary
                    nseg=length(xia);
                    xs=[xia(end),x( iis(end)+1:end ),x(1:iis(1) ),xia(1)];
                    ys=[yia(end),y( iis(end)+1:end ),y(1:iis(1) ),yia(1)];

                    dx=abs(xs(2:end)-xs(1:end-1) +i*[ys(2:end)-ys(1:end-1)] );
                    mdx(k)=max(dx);

                    if max(dx)<DXM,
                        n0=length(pslg.x);
                        pslg.x=[pslg.x,xs];
                        pslg.y=[pslg.y,ys];
                        n1=length(pslg.x);
                        nedges=[[n0+1:n1-1];[n0+2:n1]]';
                        pslg.edges=[pslg.edges;nedges];
                        nc=nc+1;pslg.chains(nc).nodes=[n0+1:n1];
                        pslg.chains(nc).index=k;
                        pslg.chains(nc).type='starts inside, first seg';
                        pslg.chains(nc).BI=1;
                    end
                    for j=2:2:nseg-1,
                        xs=[xia(j),x( iis(j):iis(j+1) ),xia(j+1)];
                        ys=[yia(j),y( iis(j):iis(j+1) ),yia(j+1)];
                        
                        dx=abs(xs(2:end)-xs(1:end-1) +i*[ys(2:end)-ys(1:end-1)] );
                        if max(dx)<DXM,
                            n0=length(pslg.x);
                            pslg.x=[pslg.x,xs];
                            pslg.y=[pslg.y,ys];
                            n1=length(pslg.x);
                            nedges=[[n0+1:n1-1];[n0+2:n1]]';
                            pslg.edges=[pslg.edges;nedges];
                            nc=nc+1;pslg.chains(nc).nodes=[n0+1:n1];
                            pslg.chains(nc).type='starts inside';
                            pslg.chains(nc).index=k;
                            pslg.chains(nc).BI=1;
                        end
                    end
                end
                
            else %isempty(xi)--> island entirely inside bounding box
                xs=[x(1:end-1)];
                ys=[y(1:end-1)];
                dx=abs(xs(2:end)-xs(1:end-1) +i*[ys(2:end)-ys(1:end-1)] );
                if and(length(xs)>2,max(dx)<DXM) %no degenerate islands
                    n0=length(pslg.x);
                    pslg.x=[pslg.x,xs];
                    pslg.y=[pslg.y,ys];
                    n1=length(pslg.x);
                    nedges=[[n0+1:n1-1];[n0+2:n1]]';
                    nedges=[nedges;[n1,n0+1]];
                    pslg.edges=[pslg.edges;nedges];
                    nc=nc+1;pslg.chains(nc).nodes=[n0+1:n1,n0+1];
                    pslg.chains(nc).type='interior island';
                    pslg.chains(nc).index=k;
                    pslg.chains(nc).BI=0;
                end
        end
        end
    end % if n>3, ar>
     if mod(k,1000)==0,disp(['Land segments compleate: ',num2str(k/N)]);,end
end

%OK - now revisit outer boundary!
eval(['save -v7.3 ',FileOutMatlab,' pslg Blon Blat isplot']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%restart script from here if needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear
%close all
eval(['load ',FileOutMatlab]);
%load pslgOSM.mat
%%Bx=[Blon(1),Blon(2),Blon(2),Blon(1),Blon(1)]
%%By=[Blat(1),Blat(1),Blat(2),Blat(2),Blat(1)]

%make outer bound
n=length(pslg.x);
%pslg.x=[pslg.x,Bx(3:4)];%add northwest and north east box corner nodes
%pslg.y=[pslg.y,By(3:4)];
pslg.x=[pslg.x,Bx(2:4)];%add northwest and north east box corner nodes
pslg.y=[pslg.y,By(2:4)];

nc=length(pslg.chains);%add northwest and north east box corner chains
pslg.chains(nc+1).nodes=n+1;
nc=length(pslg.chains);
pslg.chains(nc+1).nodes=n+2;
nc=length(pslg.chains);
pslg.chains(nc+1).nodes=n+3;
nc=length(pslg.chains);

xx =153.0400 %south east corner of Austrilia
yy = Blat(1) %check that these are the endpoints described for j0, j1!!
[m0,j0]=min(abs(pslg.x+i*pslg.y-xx-i*yy))
if isplot,plot(xx,yy,'ro',pslg.x(j0),pslg.y(j0),'bx');end

xx =  288.3143%check that these are the endpoints described for j0, j1!!
yy =  Blat(1) %south west corner of S. America
[m1,j1]=min(abs(pslg.x+i*pslg.y-xx-i*yy))
if isplot,plot(xx,yy,'ro',pslg.x(j1),pslg.y(j1),'bx');end

nc=length(pslg.chains);
pslg.chains(nc+1).nodes=[j0,j1];
nc=length(pslg.chains);

for k=1:nc
    if mod(k,100)==0,disp(['Labeling chains compleate: ',num2str(k/nc)]);end
    spx(k)=pslg.x(pslg.chains(k).nodes(1));
    spy(k)=pslg.y(pslg.chains(k).nodes(1));
    epx(k)=pslg.x(pslg.chains(k).nodes(end));
    epy(k)=pslg.y(pslg.chains(k).nodes(end));
 end
 
%make boundary order index along boundary for start points
SN=10.*[max(abs(pslg.x))+max(abs(pslg.y))];%large number to seperate edges
c=0*spy;
Deps1=10^-10
j=find(abs(spy-By(1))<Deps1);
c(j)=SN+spx(j);
j=find(abs(spx-Bx(2))<Deps1);
c(j)=2*SN+spy(j);
j=find(abs(spy-By(3))<Deps1);
c(j)=3*SN-spx(j);
j=find(abs(spx-Bx(4))<Deps1);
c(j)=4*SN-spy(j);

%make boundary order index along boundary for end points
d=0*spy;
j=find(abs(epy-By(1))<Deps1);
d(j)=1*SN+epx(j);
j=find(abs(epx-Bx(2))<Deps1);
d(j)=2*SN+epy(j);
j=find(abs(epy-By(3))<Deps1);
d(j)=3*SN-epx(j);
j=find(abs(epx-Bx(4))<Deps1);
d(j)=4*SN-epy(j);


if isplot==1,figure;plot(pslg.x,pslg.y,'k.');hold on;end
epi=j0;%start at lower left boundary corner
obc=[nc];%start at lower left boundary corner
v=SN+pslg.x(epi);
j=find(c>v);
[mm,m]=min(c(j));
l=j(m);
epi=[epi,pslg.chains(l).nodes];
obc=[obc,l];
n=0;
%while(epi(end)~=epi(1))
while isempty( find(epi(2:end)==epi(1)) )
    v=d(l);
    j=find(c>v);
    [mm,m]=min(c(j));
    l=j(m)
    epi=[epi,pslg.chains(l).nodes];
    obc=[obc,l];
  if isplot==1,  plot(pslg.x(epi),pslg.y(epi),'r.-');pause(.001);end
end
j=find(epi(2:end)==epi(1));
epi=epi(1:(j+1));
%vvvvvvvvvvv Add fixed grid points on open boundary
% to prevent "curved" boundaries due to projection details
xob=pslg.x(epi);
yob=pslg.y(epi);
dob=abs([xob(2:end)-xob(1:end-1)]+i*[yob(2:end)-yob(1:end-1)]);
dob=[dob,abs([xob(end)-xob(1)]+i*[yob(1)-yob(end)])];
nb=length(xob);
dmin=1;%node spacing around boundary
epit=epi;
for k=1:nb-1
    d=abs([xob(k+1)-xob(k)]+i*[yob(k+1)-yob(k)]);
    if d>dmin
        npl=round(d/dmin)-1;
        dx=xob(k+1)-xob(k);
        dy=yob(k+1)-yob(k);
        xp=xob(k)+dx*[1:npl-1]/npl;
        yp=yob(k)+dy*[1:npl-1]/npl;
        np=length(xp);
        n=length(pslg.x);
        pslg.x=[pslg.x,xp];
        pslg.y=[pslg.y,yp];
        NewEdges=[[epi(k),n+1];[[n+1:n+np-1]',[n+2:n+np]'];[n+np,epi(k+1)]];
        %pslg.edges=[pslg.edges;NewEdges];
        j=find(epi(k)==epit );
        epit=[epit(1:j),n+1:n+np,epit(j+1:end)];
    end
end
epi=epit;
% to prevent "curved" boundaries due to projection details
%^^^^^^^^^^^^ Add fixed grid points on open boundary


% add edges in obc
epiv=epi(:);
NOE=length(epiv)
OutterEdges=[epiv(1:NOE-1),epiv(2:NOE)];
pslg.edges=[pslg.edges;OutterEdges];

eval(['save -v7.3 ',FileOutMatlab,' pslg Blon Blat isplot']);

nodelist=[];
for k=1:nc
    if mod(k,100)==0,disp(['Finding interior chains compleate: ',num2str(k/nc)]);,end
    if pslg.chains(k).nodes(1)==pslg.chains(k).nodes(end)
        n=pslg.chains(k).nodes(1);
        [inpoly,onpoly]=insidepoly( pslg.x(n), pslg.y(n),xob,yob);
        if and(inpoly==1,onpoly==0)
           nodelist=union(nodelist,pslg.chains(k).nodes);
        end
    end
end

eval(['save -v7.3 ',FileOutMatlab,' pslg Blon Blat isplot nodelist epi']);

nodelistXB=union(nodelist,epi);
pslgb=subpslgFast(pslg,nodelistXB);

%remove "random" duplicate nodes
%z=pslgb.x+i*pslgb.y;
%[zu,j,k]=unique(z);
%pslgc=subpslgFast(pslgb,j);
pslgc=pslgb

%remove duplicate edges
pslgc.edgesS=sort(pslgc.edges')';
pslgc.edgesSU=unique(pslgc.edgesS,'rows')
pslgc.edges=pslgc.edgesSU;

pslg=pslgc

eval(['save -v7.3 ',FileOutMatlab,' pslg Blon Blat isplot']);

%save PSLG to jigsaw .msh format
geom=pslg2geom(pslg)
savemsh(FileOutJigsaw,geom)

