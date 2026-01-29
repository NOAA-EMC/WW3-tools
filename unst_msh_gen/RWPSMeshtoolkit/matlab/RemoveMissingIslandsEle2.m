function gnew=RemoveMissingIslandsEle2(g,ax,p,MinBndDist,isplot);
%function gnew=HandEditMesh(g,ax,p);
% Remove parts of mesh, g inside closed curves defined in pslg
% p.  The action is preformed only within the axis ax. The
% frame work is intended to expand to moving nodes, removing elements, etc
%This wont touch nodes within MinBndDist (m) of existing boundary
if nargin<4
    MinBndDist=1000;
end
if nargin <5,
    isplot=0;
end

xc=mean(g.x(g.e'));
yc=mean(g.y(g.e'));

jx=find(and(xc>ax(1),xc<ax(2)));
jy=find(and(yc>ax(3),yc<ax(4)));
jg=intersect(jx,jy);

jx=find(and(p.x>ax(1),p.x<ax(2)));
jy=find(and(p.y>ax(3),p.y<ax(4)));
jp=intersect(jx,jy);
p0=subpslgFast(p,jp);
chains=edges2chains(p0.edges);

nc=length(chains);
jb0=[];
for k=1:nc
    k/nc
    n=chains(k).nodes;
    if length(n)>2,
        if n(1)==n(end),
            jins=find(inside(xc(jg),yc(jg),p0.x(n),p0.y(n)));
            nin=g.e(jg(jins),:);
            X=g.x(nin);
            Y=g.y(nin);
            jnins=find( sum(inside(X,Y,p0.x(n),p0.y(n) )')==3);
            jb0=[jb0;jnins(:)];
        end
    end
end





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
if isplot,
    clf;ph=patch(x(e'),y(e'),z(e'));cm=colormap('jet');shading interp;colorbar;axis equal;
    hold on;plot(p0.x(p0.edges'),p0.y(p0.edges'),'k.-');
    set(ph,'EdgeColor','w');set(ph,'EdgeAlpha',.2);
end

nc=length(chains);
jb0=[];
for k=1:nc
    n=chains(k).nodes;
    if length(n)>2,
        if n(1)==n(end),
            jins=find(inside(x,y,p0.x(n),p0.y(n)));
            jb0=[jb0;jins(:)];
        end
    end
end

lat2m=110574.;
lon2m=111320.;

jb0X=[];
for k=1:length(jb0);
    n=jb0(k);

    d=min(abs( lon2m*( x(n)-xb )*cos(pi*y(n)/180) + ...
             i*lat2m*( y(n)-yb ) ));
    if d>MinBndDist,jb0X=[jb0X,n];end
end
jb0=jb0X;

if isplot,
    hold on
    plot(x(jb0),y(jb0),'ro');
end

[g1,kr0]=submeshFastEle(g0,jb0);

if isplot,
    figure;
    x=g1.x;y=g1.y;z=g1.z;e=g1.e;
    clf;ph=patch(x(e'),y(e'),z(e'));cm=colormap('jet');shading interp;colorbar;axis equal;
    hold on;plot(p0.x(p0.edges'),p0.y(p0.edges'),'k.-');
    set(ph,'EdgeColor','w');set(ph,'EdgeAlpha',.2);
end
%jb=jg(jb0);
jb=jg(kr0);
nn=length(g.x);
gnew=submeshFast(g,setdiff(1:nn,jb)); %only removes nodes 


