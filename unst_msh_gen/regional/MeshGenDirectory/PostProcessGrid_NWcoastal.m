

SetPath

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Input file from jigsaw and boundary file used in it's creation:

SetPath
isplot=0;
outdir='NWcoastal/'
pslgfile=PSLGfile, %global variable filename set in SetPath

jigsawout='RWPS.F.LLH'

WW3FileOut='NWcoastal.WW3.msh'

%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP (A) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=loadmsh([outdir,jigsawout,'.msh']);

%remove sand points on boundary

RemoveSandPoints([outdir,jigsawout,'.msh'],pslgfile,[outdir,jigsawout,'.NSP.msh'],[outdir,jigsawout,'.NSP.WW3.msh']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP (B) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g=loadmshWW3( [outdir,jigsawout,'.NSP.WW3.msh']);
g.x=LonCon(g.x);

if isplot,
    LS=ComputeLengthScale_wgs84_MEL(x,y,e);LSn=Ele2Nodes(x,y,e,LS);
    clf;ph=patch(x(e'),y(e'),LSn(e'));cm=colormap('jet');shading interp;axis equal;
    caxis([0,12])
    colormap(flip(cm));
end

WriteWW3MeshX(g,'RWPS.WW3a.msh');


%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP (C) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%clear
g=loadmshWW3('RWPS.WW3a.msh')
%Remove key islands remaining in mesh that are present in the PSLG file

p=loadmsh(pslgfile);

p.x=p.point.coord(:,1);
p.y=p.point.coord(:,2);
p.edges=p.edge2.index(:,1:2);
p.x=LonCon(p.x);
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
dxb=5
%close all
RemoveIsl=1;
n=1
if RemoveIsl,
   
    for xx=(ceil(min(g.x))+dxb):dx:(floor(max(g.x))-dxb)
        ax11=[xx-dx,xx+dx,min(g.y)+dxb,max(g.y)-dxb]
        N(n)=length(gnew.x);
%	ax11=[min(g.x)-dxb,max(g.x)+dxb, min(g.y)-dxb,max(g.y)+dxb]
        AX(n,:)=ax11;
        gnew=RemoveMissingIslandsEle(gnew,ax11,p,MinBndDist,0);
        n=n+1;
        length(gnew.x)
        %figure(3);clf;plot(N,'ko-');pause(.001)
    end

end

g=gnew

%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP (D) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g=RemoveSandPointsWW3(g,'RWPS.WW3a1.msh')

%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP (H) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fix boundary warping from projection in mesh generation
%clear
close all
g=loadmshWW3('RWPS.WW3a1.msh')
g0=g;

Blon=[round(min(g.x)),round(max(g.x))]
Blat=[round(min(g.y)),round(max(g.y))]
Blon=LonCon(Blon);

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
dx=.005;
js=find(and( g.x(jb)>min(Blon),g.x(jb)<min(Blon)+dx  ));
g.x(jb(js))=min(Blon);
if isplot,plot(g.x(jb(js)),g.y(jb(js)),'r.');end
js=find(and( g.x(jb)<max(Blon),g.x(jb)>max(Blon)-dx  ));
g.x(jb(js))=max(Blon);
if isplot,plot(g.x(jb(js)),g.y(jb(js)),'r.');end
WriteWW3MeshX(g,'RWPS.WW3b.msh')

%confirm no introduction of sand points
%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP (I) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g0=RemoveSandPointsWW3(g,'RWPS.WW3c.msh')

%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP (J) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WriteWW3MeshX(g0,WW3FileOut);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Make plots of of the mesh and some input fields
MakeFiguresNWCoastal
