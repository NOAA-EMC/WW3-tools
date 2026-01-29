function ConnectDisjointPSLGs(flin1,flin2,flout)
%function ConnectDisjointPSLGs(flin1,flin2,flout)
% Join 2 PSLGs together with shortest lines.
% flin1= PSLG saved in jigsaw "geom" format (big pslg)
% flin2= PSLG saved in jigsaw "geom" format (little(pslg))
% flout= where to save joined PSLG
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
