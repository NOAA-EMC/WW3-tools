function gx=ConnectDisjointMesh(g,g0,flout,N)
%function ConnectDisjointPSLGs(flin1,flin2,flout)
% Join 2 PSLGs together with shortest lines.
% flin1= PSLG saved in jigsaw "geom" format (big pslg)
% flin2= PSLG saved in jigsaw "geom" format (little(pslg))
% flout= where to save joined PSLG
Smin=500;

g=loadmsh('RWPSMeshFilesPIXNO/RWPS.F.LLH.msh')
g.x=g.point.coord(:,1);g.y=g.point.coord(:,2);g.z=g.point.coord(:,3);g.e=g.tria3.index(:,1:3);
g0=loadmsh('../RWPSLakes/RWPS.Sebago.msh')
g0.x=g0.point.coord(:,1);g0.y=g0.point.coord(:,2);g0.z=g0.point.coord(:,3);g0.e=g0.tria3.index(:,1:3);


bnd=detbndy(g.e);
bndu=unique(bnd);

xb=g.x(bndu);
yb=g.y(bndu);

bnd0=detbndy(g0.e);
bndu0=unique(bnd0);
xb0=g0.x(bndu0);
yb0=g0.y(bndu0);


lon2m=111320.
lat2m=110574.
clear D J
for k=1:length(bndu0)
    [D(k),J(k)]=min(abs(  lon2m*cos(pi*yb0(k)/180).*(xb0(k)-xb) + i*lat2m*[yb0(k)-yb] ) );
end
[D0,k0]=min(D);%closest node lake
k=J(k0);%closest node global

Nseg=round(D0/Smin)
n0=bndu0(k0)
n=bndu(k)

j0=find(bnd0(:,1)==n0);
m0=bnd0(j0,2);
j=find(bnd(:,1)==n);
m=bnd(j,2);

xng=g.x(n);
yng=g.y(n);
xng=g.x(n);
yng=g.y(n);
xp=g.x(n)+(g0.x(n0)-g.x(n))*(1:Nseg-1)/(Nseg);
yp=g.y(n)+(g0.y(n0)-g.y(n))*(1:Nseg-1)/(Nseg);

xq=g.x(m)+(g0.x(m0)-g.x(m))*(1:Nseg-1)/(Nseg);
yq=g.y(m)+(g0.y(m0)-g.y(m))*(1:Nseg-1)/(Nseg);
clf;plot(xb,yb,'k.',xb0,yb0,'b.',xp,yp,'r.',xq,yq,'c.')

g1.x=[g.x(:);g0.x(:)];
g1.y=[g.y(:);g0.y(:)];
g1.e=[g.e;length(g.x)+g0.e];

Nj=length(g1.x);
[NE,three]=size(g1.e);
g1.x=[g1.x(:);xp];
g1.y=[g1.y(:);yp];
Np=length(g1.x);
g1.x=[g1.x(:);xq];
g1.y=[g1.y(:);yq];
Nq=length(g1.x);





[ils,jls]=find(ps.edges==ks);
ksN=ps.edges(ils,:);
ksNu=unique(ksN);
ksNux=setdiff(ksNu,ks);
[ne,two]=size(ps.edges);
ps.edges=ps.edges(setdiff(1:ne,ils),:);%remove edges to close node

[ilg,jlg]=find(pg.edges==kg);
kgN=pg.edges(ilg,:);
kgNu=unique(kgN);
kgNux=setdiff(kgNu,kg);
[ne,two]=size(pg.edges);
pg.edges=pg.edges(setdiff(1:ne,ilg),:);%remove edges to close node


for k=1:length(bndu0)
    D(k)=min(abs)


pg=loadmsh(flin1)
pg.x=pg.point.coord(:,1);
pg.y=pg.point.coord(:,2);
pg.edges=pg.edge2.index(:,1:2);

ps=loadmsh(flin2);
ps.x=ps.point.coord(:,1)+360;
ps.y=ps.point.coord(:,2);
ps.edges=ps.edge2.index(:,1:2);

lon2m=111320.
lat2m=110574.
clear D
for k=1:length(ps.x)
    D(k)=min(abs(  lon2m*cos(pi*ps.y(k)/180).*(ps.x(k)-pg.x) + i*lat2m*[ps.y(k)-pg.y] ) );
end
[ms,ks]=min(D);
[mg,kg]=min(abs(  lon2m*cos(pi*ps.y(ks)/180).*(ps.x(ks)-pg.x) + i*lat2m*[ps.y(ks)-pg.y] ) );
%find neighbors of ks

[ils,jls]=find(ps.edges==ks);
ksN=ps.edges(ils,:);
ksNu=unique(ksN);
ksNux=setdiff(ksNu,ks);
[ne,two]=size(ps.edges);
ps.edges=ps.edges(setdiff(1:ne,ils),:);%remove edges to close node

[ilg,jlg]=find(pg.edges==kg);
kgN=pg.edges(ilg,:);
kgNu=unique(kgN);
kgNux=setdiff(kgNu,kg);
[ne,two]=size(pg.edges);
pg.edges=pg.edges(setdiff(1:ne,ilg),:);%remove edges to close node

nn=length(pg.x);
p=joinpslg(pg,ps);
[ne,two]=size(p.edges);
p.edges(ne+1,:)=[nn+ksNux(1),kgNux(1)];
p.edges(ne+2,:)=[nn+ksNux(2),kgNux(2)];

clf
plot(pg.x,pg.y,'k.',ps.x,ps.y,'b.')
hold on;
je=ne+1:ne+2;
plot(p.x(p.edges(je,:)'), p.y(p.edges(je,:)'),'go-');

geom=pslg2geom(p)
savemsh(flout,geom)
