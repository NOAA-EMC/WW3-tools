function OpenBndNodes=FindOuterBndNoPSLG(meshflin,ComputeBndy)
% function OpenBndNodes=FindOuterBndNoPSLG(meshflin,ComputeBndy)
% find the open boundary nodes from jigsaw format .msh file
% representing an unstructured mesh, and find the estimated 
% open ocean boundary nodes assuming the mesh is on an oriented 
% rectangle.
%
%   inputs:
%       meshflin : jigsaw format .msh file representing an unstructured mesh
%       ComputeBndy : ComputeBndy =1 in all cases
%   output: 
%       OpenBndNodes : list of open ocean boundary nodes numbers
%

Dmin=1000; %critical distance in meters to find boundary nodes

g=loadmsh(meshflin)

e=g.tria3.index(:,1:3);
x=g.point.coord(:,1);y=g.point.coord(:,2);f=g.point.coord(:,3);

%assume rectangular mesh 
x0=min(x);x1=max(x);
y0=min(y);y1=max(y);

lat2m=110574.
lon2m=111320.
    
DIE=[lon2m*(x-x0).*cos(pi*y/180),...
    lon2m*(x1-x).*cos(pi*y/180),...
    lat2m*(y-y0),lat2m*(y1-y)];
DIE=abs(DIE);
D=min(DIE');

jc=find(D<Dmin);%within 100m of edge
OpenBndNodes=jc;

%if the above isn't discriminating islands try:
if nargin>1
    if(ComputeBndy)
        bnd=detbndy(e);%find boundary elements including islands
        jb=unique(bnd(:))
        OpenBndNodes=intersect(jb,jc);
    end
end
