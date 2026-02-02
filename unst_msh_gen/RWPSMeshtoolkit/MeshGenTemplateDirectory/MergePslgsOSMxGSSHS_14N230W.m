

OSMmsh='../RWPSmeshtrial.OSM/PSLGboundary1kmP_NewOrleans.msh'
GSHHGmsh='../RWPSmeshtrial.PIXX.AP1km/PSLGboundary1kmP_NewOrleans.msh'
p0=loadmsh(OSMmsh)
p1=loadmsh(GSHHGmsh)

p0.x=p0.point.coord(:,1);
p0.y=p0.point.coord(:,2);
p0.edges=p0.edge2.index(:,1:2);


p1.x=p1.point.coord(:,1);
p1.y=p1.point.coord(:,2);
p1.edges=p1.edge2.index(:,1:2);

close all
plot(p0.x,p0.y,'b.',p1.x,p1.y,'r.')

yc=14.5 
xc=230
jx0=find(p0.x<xc);
jy0=find(p0.y<yc);
j0=intersect(jx0,jy0);

jx1=find(p1.x<xc);
jy1=find(p1.y<yc);
j1=intersect(jx1,jy1);

p0a=subpslgFast(p0,j0);

n1=length(p1.x);
k1=setdiff(1:n1,j1);
p1a=subpslgFast(p1,k1);

close all
plot(p0a.x,p0a.y,'b.',p1a.x,p1a.y,'r.')
hold on
plot([xc,xc],[min(p1.y),max(p1.y)],'c')
plot([min(p1.x),max(p1.x)],[yc,yc],'c')


p=joinpslg(p0a,p1a)


j=find(p.x<min(p.x)+.01);%Western Bounadry
na=find(p.y(j)>yc  )
nb=find(p.y(j)<yc  )
[yabove,ia]=min( p.y(j(na) ))
[ybelow,ib]=max( p.y(j(nb) ))
ma=j(na(ia))
mb=j(nb(ib))

plot(p.x(ma),p.y(ma),'gx',p.x(mb),p.y(mb),'mx')

p.edges=[p.edges;[ma,mb]];%knit together West Boundary


j=find(p.y<min(p.y)+.01);%Southern Bounadry
na=find(p.x(j)>xc  )
nb=find(p.x(j)<xc  )
[ywest,ia]=min( p.x(j(na) ))
[yeast,ib]=max( p.x(j(nb) ))
ma=j(na(ia))
mb=j(nb(ib))
plot(p.x(ma),p.y(ma),'gx',p.x(mb),p.y(mb),'mx')
%p.edges=[p.edges;[ma,mb]];%knit together South Boundary
p.edges=[p.edges;[mb,ma]];%knit together South Boundary

pslg=p;
save pslgOSMxGSHHS.mat pslg

geom=pslg2geom(pslg)

savemsh('PSLGboundaryOSMxGSHHS1km.msh',geom)


