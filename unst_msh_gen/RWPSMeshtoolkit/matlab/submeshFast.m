function h=submeshFast(g,jGoodNodes)

jGoodNodes=jGoodNodes(:)';%row format

nn=length(g.x)
[ne,three]=size(g.e)

jGoodNodes=sort(jGoodNodes);
h.x=g.x(jGoodNodes);
h.y=g.y(jGoodNodes);
h.z=g.z(jGoodNodes);

A=ismember(g.e,jGoodNodes);
k=find(sum(A')==3)';
e=g.e(k,:);
[kk,M]=ismember(1:nn,jGoodNodes);
h.e=M(e);
