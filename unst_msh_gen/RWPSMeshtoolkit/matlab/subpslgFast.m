function h=subpslgFast(g,jGoodNodes)
%function h=subpslgFast(g,jGoodNodes)
%
%Make the subset Planer Staight Line Graph(PSLG) consisting only of nodes 
% with indexes in jGoodNodes (only edges consisting spanning nodes in 
% jGoodNodes are preserved).  
%
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%Keston Smith 2022
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

jGoodNodes=jGoodNodes(:)';%row format

nn=length(g.x);
[ne,two]=size(g.edges);

jGoodNodes=sort(jGoodNodes);
h.x=g.x(jGoodNodes);
h.y=g.y(jGoodNodes);

A=ismember(g.edges,jGoodNodes);
k=find(sum(A')==2)';
edges=g.edges(k,:);
[kk,M]=ismember(1:nn,jGoodNodes);
h.edges=M(edges);
