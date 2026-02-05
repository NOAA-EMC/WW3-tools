function [h,ki]=remove_dead_nodes(g)
%function [h,ki]=remove_dead_nodes(g)
%
% remove nodes from mesh structure g that are not present in any elements
% the new mesh structure is returned as h
%
%          g.x : longitute
%          g.y : latitude
%          g.z : bathymetric depth 
%          g.e : (ne x 3) element list
%

ju=sort(unique(g.e(:)));
k=1:length(g.x);
ki=sort(setdiff(k,ju));
h.x=g.x(ju);
h.y=g.y(ju);
h.z=g.z(ju);

h.e=g.e;

for k=1:length(ki)
    j=find(g.e > ki(k) );
    h.e(j)=h.e(j)-1;
end
