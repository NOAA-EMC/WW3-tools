function [h,k]=submeshFastEle(g,jBadNodes)
%
% Make the subset of unstructured mesh g consisting only of
% elements with no nodes in jBadNodes 
%
%   input:
%          g : FE mesh structure with fields
%               g.x : longitute
%               g.y : latitude
%               g.z : bathymetric depth 
%               g.e : (ne x 3) element list
%          jBadNodes : list of nodes to throw out
%
%   output:
%          h : portion of FE mesh structure with no elements made entirely of
%              nodes in list jBadNodes
%

[ne,three]=size(g.e);
A=ismember(g.e,jBadNodes);
k=find(sum(A')==3)';%elements with nodes lying entirely within feature

elist=setdiff(1:ne,k);
g.e=g.e(elist,:);%remove elements entirely within feature;
[h,k]=remove_dead_nodes(g);%remove nodes no longer in any elements
