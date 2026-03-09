function h=subpslgFast(g,jGoodNodes)
%function h=subpslgFast(g,jGoodNodes)
%
%Make the subset Planer Staight Line Graph(PSLG) consisting only of nodes 
% with indexes in jGoodNodes (only edges consisting spanning nodes in 
% jGoodNodes are preserved).  
%   inputs:
%          g: a Planar Straight Line Graph (pslg) structure with fields
%               g.x : (n x 1) x coordinates of points in pslg
%               g.y : (n x 1 )y coordinates of points in pslg
%               g.edges : (nedges x 2) list of edges between points
%          jGoodNodes: list of points to preserve in out put h (only edges 
%                      consisting points in jGoodNodes are preserved in h
%   outputs:
%          h: subset PSLG of g consisiting only of points jGoodNodes and edges
%             spanning points in jGoodNodes
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
