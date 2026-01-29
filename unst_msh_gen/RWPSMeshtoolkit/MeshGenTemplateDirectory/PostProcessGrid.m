
outdir='RWPSMeshOSMxGSHHS.BoxesFiles/'
pslgfile='PSLGboundaryOSMxGSHHS1kmBOXES.msh'
jigsawout='RWPS.F.LLH'

%g=loadmsh('RWPSMeshOSMxGSHHS.BoxesFiles/RWPS.F.LLH.msh')
g=loadmsh([outdir,jigsawout,'.msh']);

%remove sand points on boundary
%RemoveSandPoints('RWPSMeshOSMxGSHHS.BoxesFiles/RWPS.F.LLH.msh','PSLGboundaryOSMxGSHHS1kmBOXES.msh',...
%    'RWPSMeshOSMxGSHHS.BoxesFiles/RWPS.F.LLH.NSP.msh','RWPSMeshOSMxGSHHS.BoxesFiles/RWPS.F.LLH.NSP.WW3.msh');%ileOutJigsawMesh,FileOutWW3)
RemoveSandPoints([outdir,jigsawout,'.msh'],pslgfile,[outdir,jigsawout,'.NSP.msh'],[outdir,jigsawout,'.NSP.WW3.msh']);

%Add Lakes to Mesh
%g=loadmshWW3('RWPSMeshOSMxGSHHS.BoxesFiles/RWPS.F.LLH.NSP.WW3.msh');
g=loadmshWW3( [outdir,jigsawout,'.NSP.WW3.msh']);

gS=loadmshWW3('../RWPSLakes/Sebago.NWPS.WW3.msh');
gS.x=gS.x-360;
gW=loadmshWW3('../RWPSLakes/Winnipesaukee.NWPS.WW3.msh')
gW.x=gW.x-360;
gO=loadmshWW3('../RWPSLakes/Okeechobee.NWPS.WW3.msh')
gO.x=gO.x-360;

g=CombineMesh(g,gO);
g=CombineMesh(g,gS);
g=CombineMesh(g,gW);

x=g.x;y=g.y;z=g.z;e=g.e;
LS=ComputeLengthScale_wgs84_MEL(x,y,e);LSn=Ele2Nodes(x,y,e,LS);

clf;ph=patch(x(e'),y(e'),LSn(e'));cm=colormap('jet');shading interp;axis equal;
caxis([0,12])
colormap(flip(cm));

WriteWW3MeshX(g,'RWPS.WW3a.lakes.msh');


%clear
g=loadmshWW3('RWPS.WW3a.lakes.msh')
%Remove key islands remaining in mesh
%load PSLGboundary1kmP_NewOrleans.mat
%p=loadmsh('PSLGboundaryOSMxGSHHS1kmBOXES.msh')
p=loadmsh(pslgfile);

p.x=p.point.coord(:,1);
p.y=p.point.coord(:,2);
p.edges=p.edge2.index(:,1:2);
p.x=p.x-360;
x=g.x;y=g.y;z=g.z;e=g.e;
clf;ph=patch(x(e'),y(e'),z(e'));cm=colormap('jet');shading interp;axis equal;
hold on
plot(p.x,p.y,'k.')

MinBndDist=1000;
gnew=g;

%go through all longitude and remove nodes 
dx=1;
%close all
n=1;
if 1,
    for xx=(ceil(min(g.x))+2*dx):dx:(floor(max(g.x))-2*dx)
        ax11=[xx-dx,xx+dx,min(g.y)+dx,max(g.y)-dx]
        N(n)=length(gnew.x);
        gnew=RemoveMissingIslandsEle(gnew,ax11,p,MinBndDist,0);
        n=n+1;
        n/243
        length(gnew.x)
        figure(3);clf;plot(N,'ko-');pause(.001)
    end
end
    

g=gnew
save NoPacIslAll.mat g pslgfile p




clear
load NoPacIslAll.mat

x=g.x;y=g.y;z=g.z;e=g.e;
LS=ComputeLengthScale_wgs84_MEL(x,y,e);LSn=Ele2Nodes(x,y,e,LS);

clf;ph=patch(x(e'),y(e'),LSn(e'));cm=colormap('jet');shading interp;axis equal;
caxis([0,12])
colormap(flip(cm));
%p=loadmsh('PSLGboundaryOSMxGSHHS1km.msh');
p=loadmsh(pslgfile);
p.x=p.point.coord(:,1);p.x=p.x-360;p.y=p.point.coord(:,2);p.edges=p.edge2.index(:,1:2);
hold on;
plot(p.x,p.y,'k.');
g0=g;
gnew=g
axAS =[ -179.6738 -159.5513  -18.4086   -7.1731]
%gnew=RemoveMissingIslandsEle(gnew,axAS,p,MinBndDist);
axis(axAS);pause(5)

axTofol =[ -197.5575 -196.5479    4.8724    5.7902]
axis(axTofol);pause(5)

axPA = [-177.6646 -175.5927   -0.2906    1.3435]
%gnew=RemoveMissingIslandsEle(gnew,axPA,p,MinBndDist);
axis(axPA);pause(5)

axJA =[ -172.1141 -165.7656   14.1706   19.1777]
%gnew=RemoveMissingIslandsEle(gnew,axJA,p,MinBndDist);
axis(axJA);pause(5)

axMWI =[ -179.5965 -176.0019   27.0046   29.8397]
%gnew=RemoveMissingIslandsEle(gnew,axMWI,p,MinBndDist);
axis(axMWI);pause(5)

axTin=[-215.1401 -213.4623   14.3922   15.7155]
%gnew=RemoveMissingIslandsEle(gnew,axTin,p,MinBndDist);
axis(axTin);pause(5)

axNNMI =[ -215.9834 -213.8052   19.3525   21.0704];
%gnew=RemoveMissingIslandsEle(gnew,axNNMI,p,MinBndDist);
axis(axNNMI);pause(5)

axNNMI1=[-214.8924 -213.6875   14.6480   15.5983];
%gnew=RemoveMissingIslandsEle(gnew,axNNMI1,p,MinBndDist);
axis(axNNMI1);pause(5)

axMaj =[ -177.2629 -175.7651   -0.1340    1.0474]
%gnew=RemoveMissingIslandsEle(gnew,axMaj,p,MinBndDist);
axis(axMaj);pause(5)

axHI =[ -179.1073 -153.7595   11.4935   31.4855]
%gnew=RemoveMissingIslandsEle(gnew,axHI,p,MinBndDist);
axis(axHI);pause(5)

axChuuk =[ -198.1190 -196.0417    4.7967    5.8354]
%gnew=RemoveMissingIslandsEle(gnew,axChuuk,p,MinBndDist);
axis(axChuuk);pause(5)

WriteWW3MeshX(gnew,'RWPS.WW3b.lakes.msh') %Fails


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fix new orleans
clear
g=loadmshWW3('RWPS.WW3b.lakes.msh')
x=g.x;y=g.y;z=g.z;e=g.e;
LS=ComputeLengthScale_wgs84_MEL(x,y,e);
LSn=Ele2Nodes(x,y,e,LS);
clf;ph=patch(x(e'),y(e'),LSn(e'));cm=colormap('jet');shading interp;axis equal;
colormap(flip(cm));
caxis([0,12]);

S=shaperead('../NewOrleansCoast/mz03mr26_LIX.shp')
hold on
for k=1:length(S),plot(S(k).X,S(k).Y,'k');end

axNOZ =[  -91.6701  -87.7499   28.7284   30.9173]
axis(axNOZ)
gnew=g;
%if nescesary run:
display('enter 1 if you want to edit the mesh around New Orleans' )
gnew=RemoveMeshParts(g,axNOZ,S,[0,50]);

WriteWW3MeshX(gnew,'RWPS.WW3c.lakes.msh')
%confirm no introduction of sand points
g=RemoveSandPointsWW3(gnew,'RWPS.WW3d.lakes.msh')
WriteWW3MeshX(g,'RWPS.WW3d.lakes.msh')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fix boundary warping from projection in mesh generation
clear
g=loadmshWW3('RWPS.WW3d.lakes.msh')
Blon=[129.91 10.71]
Blat=[-30.42 79.99]
Blon(1)=Blon(1)-360
x=g.x;y=g.y;z=g.z;e=g.e;
jb=g.bnd;

clf;
plot(g.x(jb),g.y(jb),'k.');hold on
js=find(g.y(jb)<min(Blat));
g.y(jb(js))=min(Blat);
plot(g.x(jb(js)),g.y(jb(js)),'r.')
js=find(g.y(jb)>max(Blat));
g.y(jb(js))=max(Blat);
plot(g.x(jb(js)),g.y(jb(js)),'r.')

dx=.5;
js=find(and( g.x(jb)>min(Blon),g.x(jb)<min(Blon)+dx  ));
g.x(jb(js))=min(Blon);
plot(g.x(jb(js)),g.y(jb(js)),'r.')
js=find(and( g.x(jb)<max(Blon),g.x(jb)>max(Blon)-dx  ));
g.x(jb(js))=max(Blon);
plot(g.x(jb(js)),g.y(jb(js)),'r.')
WriteWW3MeshX(g,'RWPS.WW3e.lakes.msh')
%confirm no introduction of sand points
g0=RemoveSandPointsWW3(g,'RWPS.WW3f.lakes.msh')

WriteWW3MeshX(g0,'RWPS.WW3g.lakes.msh')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Some Plotting



clear
g=loadmshWW3('RWPS.WW3g.lakes.msh')%RWPS.PIXAllLnwps.PP.WW3c.msh');
x=g.x;y=g.y;z=g.z;e=g.e;
LS=ComputeLengthScale_wgs84_MEL(x,y,e);LSn=Ele2Nodes(x,y,e,LS);
clf;ph=patch(x(e'),y(e'),LSn(e'));cm=colormap('jet');shading interp;axis equal;
jb=g.bnd;
hold on;
plot(g.x(jb),g.y(jb),'k.');hold on


f=z;
clf;patch(x(e'),y(e'),f(e'));cm=colormap('jet');shading interp;axis equal;
colormap(flip(cm));
caxis([0,6000]);
kprint('Bathy.jpg');

clf;
    caxis([0,6000]);
colormap(flip(cm));

    hcb=colorbar('h','position',[.1,.075,.8,.025])
    %hcb.Label.String='(km)'
    set(gca,'visible','off')
kprint('BathyColorbarH.jpg');

clf;
colorbar;
colormap(flip(cm));
caxis([0,6000]);
kprint('BathyColorbarV.jpg');


LS=ComputeLengthScale_wgs84_MEL(x,y,e);
LSn=Ele2Nodes(x,y,e,LS);

clf;patch(x(e'),y(e'),LSn(e'));cm=colormap('jet');shading interp;axis equal;
colormap(flip(cm));
caxis([0,12]);
kprint('Lengthscale.jpg');


axCAR =[  -85.4065  -61.7227    8.8200   26.8661]
ax=axCAR
axis(ax);daspect([1,cos(pi*ax(3)/180),1])
kprint('Caribean.jpg')


axAK =[ -192.1450 -152.2804   47.2479   77.6231]
ax=axAK
axis(ax);daspect([1,cos(pi*ax(3)/180),1])
kprint('Alaska.jpg')

axHI =[ -179.5950 -152.9066   13.8641   34.1996]
ax=axHI
axis(ax);daspect([1,cos(pi*ax(3)/180),1])
kprint('Hawaii.jpg')


axMI =[ -220.2161 -208.6314   12.5003   21.3274];
ax=axMI
axis(ax);daspect([1,cos(pi*ax(3)/180),1])
kprint('MarianaIslands.jpg')

axAS =[ -173.3623 -168.3684  -16.4846  -12.0090]
ax=axAS
axis(ax);daspect([1,cos(pi*ax(3)/180),1])
kprint('AmericanSamoa.jpg')


axMicro =[ -210.3818 -185.8701   -3.7494   15.5832]
ax=axMicro
axis(ax);daspect([1,cos(pi*ax(3)/180),1])
kprint('Micronesia.jpg')

axWMicro =[ -226.6896 -214.3594    5.3753   15.1002]
ax=axWMicro
axis(ax);daspect([1,cos(pi*ax(3)/180),1])
kprint('WestMicronesia.jpg')


figure;
clf;
    caxis([0,12]);
colormap(flip(cm));

    hcb=colorbar('h','position',[.1,.075,.8,.025])
    hcb.Label.String='(km)'
    set(gca,'visible','off')
kprint('LengthScaleColorbarH.jpg');

clf
    caxis([0,12]);
colormap(flip(cm));

   % hcb=colorbar('v','position',[.1,.075,.8,.025])
    vcb=colorbar('v','position',[.075,.1,.025,.8])
    vcb.Label.String='(km)'
    set(gca,'visible','off')
kprint('LengthScaleColorbarV.jpg');


clf
caxis([0,6000]);
colormap(flip(cm));

    hcb=colorbar('h','position',[.1,.075,.8,.025])
    hcb.Label.String='(m)'
    set(gca,'visible','off')
kprint('BathyColorbarH.jpg');

clf
    caxis([0,6000]);
colormap(flip(cm));

   % hcb=colorbar('v','position',[.1,.075,.8,.025])
    vcb=colorbar('v','position',[.075,.1,.025,.8])
    vcb.Label.String='(m)'
    set(gca,'visible','off')
kprint('BathyColorbarV.jpg');
qi=