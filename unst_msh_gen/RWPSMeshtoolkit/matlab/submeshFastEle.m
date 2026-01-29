function [h,k]=submeshFastEle(g,jBadNodes)

[ne,three]=size(g.e);
A=ismember(g.e,jBadNodes);
k=find(sum(A')==3)';%elements with nodes lying entirely within feature

elist=setdiff(1:ne,k);
g.e=g.e(elist,:);%remove elements entirely within feature;
[h,k]=remove_dead_nodes(g);%remove nodes no longer in any elements
