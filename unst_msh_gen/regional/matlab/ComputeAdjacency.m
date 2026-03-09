function Ac=ComputeAdjacency(e)
% function Ac=ComputeAdjacency(e)
% Compute number of adjacencent elements for each node in mesh with element list e
%
% input:    e  (ne x 3) element list (ne : number of elements in mesh)
% outout:   Ac (nn x 1) number of elements adjacent to each node (nn : number of nodes in mesh)

[Ac,GR] = groupcounts(e(:));
