function gnew=RemoveMissingIslandsPoint(g,ax,p,MinBndDist);
% function gnew=RemoveMissingIslandsPoint(g,ax,p,MinBndDist);
% Not in current use. This will be removed upon further review

if nargin<4
    MinBndDist=1000;
end

jx=find(and(g.x>ax(1),g.x<ax(2)));
jy=find(and(g.y>ax(3),g.y<ax(4)));
jg=intersect(jx,jy);
g0=submeshFast(g,jg);
bnd=detbndy(g0.e);
jbnd=unique(bnd(:));

xb=g0.x(jbnd);
yb=g0.y(jbnd);

jx=find(and(p.x>ax(1),p.x<ax(2)));
jy=find(and(p.y>ax(3),p.y<ax(4)));
jp=intersect(jx,jy);
p0=subpslgFast(p,jp);

x=g0.x;y=g0.y;z=g0.z;e=g0.e;
chains=edges2chains(p0.edges);

clf;ph=patch(x(e'),y(e'),z(e'));cm=colormap('jet');shading interp;colorbar;axis equal;
hold on;plot(p0.x(p0.edges'),p0.y(p0.edges'),'k.-')
set(ph,'EdgeColor','w');set(ph,'EdgeAlpha',.2);

LS=ComputeLengthScale_wgs84_MEL(x,y,e);LSn=Ele2Nodes(x,y,e,LS);
nc=length(chains)
lat2m=110574.
lon2m=111320.
jb0=[];
for k=1:nc
    n=chains(k).nodes;
    if length(n)>2,
        dx=p.x(n(1:end-1))-p.x(n(2:end));
        dy=p.y(n(1:end-1))-p.y(n(2:end));
        d=abs(  lon2m*dx*cos(pi*mean(p.y(n))/180)+i*lat2m*(dy) );
        P=sum(d)/1000;%perimeter length of feature
        mx=mean(p0.x(n));
        my=mean(p0.y(n));
        [m,jc]=min(abs(g0.x-mx+i*[g0.y-my] ));
        P
            LSn(jc)
            
        if LSn(jc)*6 < P, %if "new island" smaller than feature- cut it
            hold on;plot(mx,my,'ro',g0.x(jc),g0.y(jc),'cx');
            jb0=[jb0;jc];
        else
            hold on;plot(mx,my,'r.',g0.x(jc),g0.y(jc),'rx');
        end
    end
end
jc
jb0
lat2m=110574.
lon2m=111320.

jb0X=[];
for k=1:length(jb0);
    n=jb0(k);

    d=min(abs( lon2m*( x(n)-xb )*cos(pi*y(n)/180) + ...
             i*lat2m*( y(n)-yb ) ));
    if d>=MinBndDist,jb0X=[jb0X,n];end
end
jb0=jb0X

hold on
plot(x(jb0),y(jb0),'ro');
xjb0=x(jb0);
yjb0=y(jb0);
nn=length(x);
jg0=setdiff(1:nn,jb0);
g1=submeshFast(g0,jg0);

figure;
x=g1.x;y=g1.y;z=g1.z;e=g1.e;
clf;ph=patch(x(e'),y(e'),z(e'));cm=colormap('jet');shading interp;colorbar;axis equal;
hold on;plot(p0.x(p0.edges'),p0.y(p0.edges'),'k.-')
set(ph,'EdgeColor','w');set(ph,'EdgeAlpha',.2);
plot(xjb0,yjb0,'r.');

jb=jg(jb0);
nn=length(g.x);
gnew=submeshFast(g,setdiff(1:nn,jb));


