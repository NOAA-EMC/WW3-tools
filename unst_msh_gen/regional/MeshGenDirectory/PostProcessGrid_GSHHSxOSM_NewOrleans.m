
%script to handle all post jigsaw mesh editing.
%Post processing steps are as follows
% (A) Remove sand points
% (B) Merge inland lakes --> 'RWPS.WW3a.lakes.msh'
% (C) Remove islands that jigsaw has meshed over 
% (D) Remove sand points that can be created (rarely) with island removal --> 'RWPS.WW3a1.lakes.msh'
% INACTIVE (E) Examine Pacific island resolution and make changes to island editing if needed --> 'RWPS.WW3b.lakes.msh'
% INACTIVE (F) Examine mesh near New Orleans where boundary has been merged to reflect new marine zones --> 'RWPS.WW3c.lakes.msh'
% INACTIVE (G) Remove sand points that can be created (rarely) with New Orleans editing --> 'RWPS.WW3d.lakes.msh'
% (H) Stretch open ocean boundary nodes back to origonal specified boundary rectangle --> 'RWPS.WW3e.lakes.msh'
%       this is an artificat of the projection used in jigsaw
% (I) Remove sand points that can be created (Should not matter) with boundary node stretching --> 'RWPS.WW3f.lakes.msh'
% (J) Write final WW3 mesh --> 'RWPS.WW3g.lakes.msh'
%       this is the mesh to run WW3 on.

SetPath
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Input file from jigsaw and boundary file used in it's creation:

isplot=0;
outdir='RWPS.GSHHSxOSM.NewOrleans/'
pslgfile=PSLGfile
jigsawout='RWPS.F.LLH'
WW3FileOut='RWPS.GSHHSxOSM.NewOrleans.WW3.msh'
ax=[-178,-154, 18,30]


%parameters for graphics plotting
%Local Plotting axis
ax(1:2)=LonCon(ax(1:2))
caxH=[0,5000],%color limits for bathymetry(m)
caxLS=[0,12] %color limits for mesh length scale (km)


%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP (A) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g=loadmsh([outdir,jigsawout,'.msh']);

PlotJigsawUnstMesh([outdir,jigsawout,'.msh'],ax,caxH,caxLS);

%remove sand points on boundary

RemoveSandPoints([outdir,jigsawout,'.msh'],pslgfile,[outdir,jigsawout,'.NSP.msh'],[outdir,jigsawout,'.NSP.WW3.msh']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP (B) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Add Lakes to Mesh
g=loadmshWW3( [outdir,jigsawout,'.NSP.WW3.msh']);

gS=loadmshWW3([LakeDir,'Sebago.NWPS.WW3.msh']);
gS.x=gS.x-360;
gW=loadmshWW3([LakeDir,'Winnipesaukee.NWPS.WW3.msh'])
gW.x=gW.x-360;
gO=loadmshWW3([LakeDir,'Okeechobee.NWPS.WW3.msh'])
gO.x=gO.x-360;

g=CombineMesh(g,gO);
g=CombineMesh(g,gS);
g=CombineMesh(g,gW);

x=g.x;y=g.y;z=g.z;e=g.e;

if isplot,
    LS=ComputeLengthScale_wgs84_MEL(x,y,e);LSn=Ele2Nodes(x,y,e,LS);
    clf;ph=patch(x(e'),y(e'),LSn(e'));cm=colormap('jet');shading interp;axis equal;
    caxis([0,12])
    colormap(flip(cm));
end

WriteWW3MeshX(g,'RWPS.WW3a.lakes.msh');


%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP (C) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%clear
g=loadmshWW3('RWPS.WW3a.lakes.msh')
%Remove key islands remaining in mesh that are present in the PSLG file

p=loadmsh(pslgfile);

p.x=p.point.coord(:,1);
p.y=p.point.coord(:,2);
p.edges=p.edge2.index(:,1:2);
p.x=p.x-360;
x=g.x;y=g.y;z=g.z;e=g.e;

if isplot,
    clf;ph=patch(x(e'),y(e'),z(e'));cm=colormap('jet');shading interp;axis equal;
    hold on
    plot(p.x,p.y,'k.')
end

MinBndDist=1000;
gnew=g;

%go through all longitude and remove nodes 
dx=2;
dxb=10
%close all
RemoveIsl=1;
n=1
if RemoveIsl,
    for xx=(ceil(min(g.x))+dxb):dx:(floor(max(g.x))-dxb)
        ax11=[xx-dx,xx+dx,min(g.y)+dxb,max(g.y)-dxb]
        N(n)=length(gnew.x);
        AX(n,:)=ax11;
        gnew=RemoveMissingIslandsEle(gnew,ax11,p,MinBndDist,0);
        n=n+1;
        n/243
        length(gnew.x)
        %figure(3);clf;plot(N,'ko-');pause(.001)
    end
%    figure(3);clf;plot(AX(:,1),N,'ko-');
end
    

g=gnew

%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP (D) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g=RemoveSandPointsWW3(g,'RWPS.WW3a1.lakes.msh')
PlotWW3Mesh('RWPS.WW3a1.lakes.msh',ax,caxH,caxLS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP (E) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p=loadmsh(pslgfile);
p.x=p.point.coord(:,1);p.x=p.x-360;p.y=p.point.coord(:,2);p.edges=p.edge2.index(:,1:2);

x=g.x;y=g.y;z=g.z;e=g.e;
if isplot,
    LS=ComputeLengthScale_wgs84_MEL(x,y,e);LSn=Ele2Nodes(x,y,e,LS);
    clf;ph=patch(x(e'),y(e'),LSn(e'));cm=colormap('jet');shading interp;axis equal;
    caxis([0,12])
    colormap(flip(cm));
    hold on;
    plot(p.x,p.y,'k.');
end

p=loadmsh(pslgfile);
p.x=p.point.coord(:,1);p.x=p.x-360;p.y=p.point.coord(:,2);p.edges=p.edge2.index(:,1:2);

g0=g;
gnew=g
if isplot
    axAS =[ -179.6738 -159.5513  -18.4086   -7.1731]
    axis(axAS);pause(5)

    axTofol =[ -197.5575 -196.5479    4.8724    5.7902]
    axis(axTofol);pause(5)

    axPA = [-177.6646 -175.5927   -0.2906    1.3435]
    axis(axPA);pause(5)

    axJA =[ -172.1141 -165.7656   14.1706   19.1777]
    axis(axJA);pause(5)

    axMWI =[ -179.5965 -176.0019   27.0046   29.8397]
    axis(axMWI);pause(5)

    axTin=[-215.1401 -213.4623   14.3922   15.7155]
    axis(axTin);pause(5)

    axNNMI =[ -215.9834 -213.8052   19.3525   21.0704];
    axis(axNNMI);pause(5)

    axNNMI1=[-214.8924 -213.6875   14.6480   15.5983];
    axis(axNNMI1);pause(5)

    axMaj =[ -177.2629 -175.7651   -0.1340    1.0474]
    axis(axMaj);pause(5)

    axHI =[ -179.1073 -153.7595   11.4935   31.4855]
    axis(axHI);pause(5)

    axChuuk =[ -198.1190 -196.0417    4.7967    5.8354]
    axis(axChuuk);pause(5)
end

WriteWW3MeshX(gnew,'RWPS.WW3b.lakes.msh')



%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP (F) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fix new orleans
%clear
g=loadmshWW3('RWPS.WW3b.lakes.msh')
%Examine and Hand Edit Mesh around new orleans to relect changes if nescesary
%Should not be nescesary with 
if isplot,
    x=g.x;y=g.y;z=g.z;e=g.e;
    LS=ComputeLengthScale_wgs84_MEL(x,y,e);
    LSn=Ele2Nodes(x,y,e,LS);
    clf;ph=patch(x(e'),y(e'),LSn(e'));cm=colormap('jet');shading interp;axis equal;
    colormap(flip(cm));
    caxis([0,12]);
    % Shape file of new marine zones in New Orleans
    S=shaperead('/scratch3/NCEPDEV/climate/Keston.Smith/RWPS/Data/us_coastline/mz03mr26_LIX.shp')
    hold on
    for k=1:length(S),plot(S(k).X,S(k).Y,'k');end

    axNOZ =[  -91.6701  -87.7499   28.7284   30.9173]
    axis(axNOZ)
    gnew=g;
    eno=input('Enter 1 if you want to edit the mesh around New Orleans:' )
    if eno,
        gnew=RemoveMeshParts(g,axNOZ,S,[0,50]);
    else
        gnew=g;
    end
else
    gnew=g;% no edit
end
   
WriteWW3MeshX(gnew,'RWPS.WW3c.lakes.msh')

%confirm no introduction of sand points
%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP (G) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=RemoveSandPointsWW3(gnew,'RWPS.WW3d.lakes.msh')


%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP (H) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fix boundary warping from projection in mesh generation
%clear
close all
g=loadmshWW3('RWPS.WW3d.lakes.msh')
g0=g;

Blon=[129.91 10.71]
Blat=[-30.42 79.99]
Blon(1)=Blon(1)-360
x=g.x;y=g.y;z=g.z;e=g.e;
jb=g.bnd;
bnd0=detbndy(g.e)

if isplot,
    clf;ph=patch(x(e'),y(e'),z(e'));
    cm=colormap('jet');shading interp;hold on
     colormap(flip(cm));
    axis equal;hold on
    plot(g.x(jb),g.y(jb),'k.');hold on
end
js=find(g.y(jb)<min(Blat));
g.y(jb(js))=min(Blat);
if isplot,plot(g.x(jb(js)),g.y(jb(js)),'r.');end
js=find(g.y(jb)>max(Blat));
g.y(jb(js))=max(Blat);
if isplot,plot(g.x(jb(js)),g.y(jb(js)),'r.');end

%dx=degree latitude delta to discriminate land points near boundary from open ocean boundary
dx=.05;
js=find(and( g.x(jb)>min(Blon),g.x(jb)<min(Blon)+dx  ));
g.x(jb(js))=min(Blon);
if isplot,plot(g.x(jb(js)),g.y(jb(js)),'r.');end
js=find(and( g.x(jb)<max(Blon),g.x(jb)>max(Blon)-dx  ));
g.x(jb(js))=max(Blon);
if isplot,plot(g.x(jb(js)),g.y(jb(js)),'r.');end
WriteWW3MeshX(g,'RWPS.WW3e.lakes.msh')

%confirm no introduction of sand points
%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP (I) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g0=RemoveSandPointsWW3(g,'RWPS.WW3f.lakes.msh')

%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP (J) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WriteWW3MeshX(g0,WW3FileOut);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Some Plotting
PlotWW3Mesh(WW3FileOut,ax,caxH,caxLS);
