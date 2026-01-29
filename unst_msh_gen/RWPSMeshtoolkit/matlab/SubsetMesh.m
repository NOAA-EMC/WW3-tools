function omesh=SubsetMesh(mesh,j)
%function omesh=SubsetMesh(mesh,j)
omesh.point.coord=mesh.point.coord(j,:);
e=mesh.tria3.index(:,1:3);

[n1,four]=size(mesh.point.coord);
[n0,four]=size(omesh.point.coord);

%S=sparse(n1,n0);
%S(1:n0,j)=1;
[A,B]=ismember(1:n1,j);
e0=B(e);%map to kept nodes and 0 for missing nodes
e0i=e0;
e0i(find(e0>0))=1;%index of nodes (1) in or (0) out
ke=find(sum(e0i')==3);
e0=e0(ke,:);

omesh.tria3.index(:,1:3)=e0;
[ne,three]=size(omesh.tria3.index);
omesh.tria3.index(:,4)=zeros(ne,1);


