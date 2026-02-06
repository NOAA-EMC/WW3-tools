function RemoveSandPoints(FileInJigsawMesh,FileInJigsawPSLG,FileOutJigsawMesh,FileOutWW3)
%function RemoveSandPoints(FileInJigsawMesh,FileInJigsawPSLG,FileOutJigsawMesh,FileOutWW3)
% Find and remove sandpoints from fininte element mesh structure g
% Because removal of sand points can create new sand points, the process is run iteratively
% within this routine. Usually, 1-3 iterations are required to clear all sandpoints in large 
% ~ 5 million node meshes generated with jigsaw. 

% A sandpoint is a boundary node with more than two adjacent boundary nodes.  Sandpoints are 
% flagged in WW3 and WW3 will not run with a mesh containing sandpoints(error in grid.out).
%
%   input:
%       FileInJigsawMesh : jigsaw format .msh file describing a finite element mesh(jigsaw output)
%       FileInJigsawPSLG : jigsaw format .msh file describing the geometry of the PSLG used 
%                          to create the mesh defined in  FileInJigsawMesh
%       FileOutJigsawMesh : jisaw format .msh file with sandpoints removed
%       FileOutWW3 : WW3 .msh file with sandpoints removed
%       

%g=loadmsh('output/RWPS.F.LLH.msh');
%g=loadmsh('../RWPSMeshTest/output.1km.10km/RWPS.F.LLH.msh');
g=loadmsh(FileInJigsawMesh)

close all;
count=0;
IsValidBoundary=0
while ~IsValidBoundary,
    e=g.tria3.index(:,1:3);
    x=g.point.coord(:,1);y=g.point.coord(:,2);f=g.point.coord(:,3);
    bnd=detbndy(e);
    bndn=unique(bnd(:));
    clear n1 n2
    for k=1:length(bndn);
        j1=find(bndn(k)==bnd(:,1));
        j2=find(bndn(k)==bnd(:,2));
        n1(k)=length(j1);
        n2(k)=length(j2);
    end
    
    nb=n1+n2;
    unique(nb)% 2 or 4

    j4=find(nb==4)
    if isempty(j4)
        IsValidBoundary=1
    end
  
    figure;
    clf;patch(x(e'),y(e'),f(e'));shading interp;colormap('jet');colorbar;
    hold on
    phsp=plot(x(bndn(j4)),y(bndn(j4)),'ro');
    title(['iteration : ',int2str(count)])
    if ~IsValidBoundary,
        count=count+1;
         bb=bndn(j4)
         bbu=unique(bb)
         display(['number of bad boundary points= ',int2str(length(bbu))])
         display(['iteration number : ',int2str(count)])
         nn=length(x);
         %jgn=setdiff([1:nn],bbu(:)');
         
         jgnBnd=setdiff([1:nn],bbu(:)');
         g.x=x(:);g.y=y(:);g.z=f(:);g.e=e;
         hBnd=submeshFast(g,jgnBnd);%can orphan nodes 
         jgnInt=unique(hBnd.e(:));%nodes still in an element
         h=submeshFast(hBnd,jgnInt);

         [ne,three]=size(h.e)
         nn=length(h.x)
         g0=g;%put back into jigsaw format
         g0=rmfield(g0,"tria3")
         g0=rmfield(g0,"point")
         ep=[h.e,zeros(ne,1)];
         g0.tria3.index=ep;
         g0.point.coord=[h.x(:),h.y(:),h.z(:),zeros(nn,1)];
         g0.value=-h.z;
         g=g0;
    end
end

savemsh(FileOutJigsawMesh,g);
WriteWW3Mesh(FileOutJigsawMesh,FileInJigsawPSLG,FileOutWW3);
