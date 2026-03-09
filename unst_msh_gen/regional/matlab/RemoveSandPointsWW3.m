function g=RemoveSandPointsWW3(g,FileOutWW3)
%function g=RemoveSandPointsWW3(g,FileOutWW3)
% Find and remove sandpoints from fininte element mesh structure g
% Because removal of sand points can create new sand points, the process is run iteratively
% within this routine. Usually, 1-3 iterations are required to clear all sandpoints in large 
% ~ 5 million node meshes generated with jigsaw. 

% A sandpoint is a boundary node with more than two adjacent boundary nodes.  Sandpoints are 
% flagged in WW3 and WW3 will not run with a mesh containing sandpoints(error in grid.out).
%
%   input:
%          g : FE mesh structure with fields
%               g.x : longitute
%               g.y : latitude
%               g.z : bathymetric depth 
%               g.e : (ne x 3) element list
%          FileOutWW : file name to write WW3 msh file for g 
%
%   outputs:
%           g : input mesh with all sandpoints removed by deletion
%       

close all;
count=0;
IsValidBoundary=0
while ~IsValidBoundary, 
    x=g.x;y=g.y;f=g.z;e=g.e;
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
    plot(x(bndn(j4)),y(bndn(j4)),'ro');
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
         g=h;
    end
end

if nargin >1,
    WriteWW3MeshX(g,FileOutWW3);
end
