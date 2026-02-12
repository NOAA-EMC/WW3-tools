% Script to merge PSLG created from GSHHS coastline with PSLG created from OSM coastline data.
% GSHHS coastline is used everywhere except with in rectangles described in rows of AX.

isplot=0
OSMmsh='./GlobalCoastlineOSM.PSLG.msh' 
GSHHGmsh='./GlobalCoastlineGSHHS.PSLG.msh'
p0=loadmsh(OSMmsh)
p1=loadmsh(GSHHGmsh)

p0.x=p0.point.coord(:,1);
p0.y=p0.point.coord(:,2);
p0.edges=p0.edge2.index(:,1:2);


p1.x=p1.point.coord(:,1);
p1.y=p1.point.coord(:,2);
p1.edges=p1.edge2.index(:,1:2);

axAS =  [183.0293  191.6089  -24.5496  -11.9295]
axPalau =[  130.8368  141.8230    3.9914   10.1256]
axSP =[  211.0219  220.1817  -18.9565  -13.8421]
axCP =[  151.1962  152.3336    6.9451    7.5801]
axGOM =[  290.8507  291.1934   43.7119   43.9032]
axCar =[  272.0041  272.7733   17.1417   17.5712]
axCarB =[  272.9356  274.2177   16.0475   16.7634];

AX=[axAS;axPalau;axSP;axCP;axGOM;axCar;axCarB]
[nax,four]=size(AX)

if isplot,
    close all
    clf;
    plot(p0.x,p0.y,'b.',p1.x,p1.y,'r.');
    axis equal;
    hold on
    for k=1:nax
        ax=AX(k,:);
        BoundingBox(ax(1:2),ax(3:4),'c');
    end
end

p=p1;%gshhs baseline
for k=1:nax
    j0=FindPointsAx(AX(k,:),p0.x,p0.y);
    j=FindPointsAx(AX(k,:),p.x,p.y);
    p0a=subpslgFast(p0,j0);
    nn=length(p.x);
    k=setdiff(1:nn,j);
    p=subpslgFast(p,k);%remove gshhs features in box
    p0a=subpslgFast(p0,j0);%get OSM features in box
    p=joinpslg(p,p0a);%add OSM features to remaining gshhs 
end

pslg=p;
save pslgOSMxGSHHS.BOXES.mat pslg

geom=pslg2geom(pslg)

savemsh('PSLGboundaryOSMxGSHHS.BOXES.msh',geom)

hold on;
plot(p.x,p.y,'k.');
title('red-gshhs, blue- OSM, black- final')
