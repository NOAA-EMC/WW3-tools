function h=submeshFast(g,jGoodNodes)
%
% Make the subset of unstructured mesh g consisting only of nodes 
% with indexes in jGoodNodes (and only elements consisting entirely
% of nodes in jGoodNodes are preserved).  
%   input:
%          g : FE mesh structure with fields
%               g.x : longitute
%               g.y : latitude
%               g.z : bathymetric depth 
%               g.e : (ne x 3) element list
%          jGoodNodes : list of nodes to keep
%
%   output:
%          h : portion of FE mesh structure,g spanned by nodes jGoodNodes
%

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
