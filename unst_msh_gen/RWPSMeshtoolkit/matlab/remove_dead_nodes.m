function [h,ki]=remove_dead_nodes(g)


ju=sort(unique(g.e(:)));
k=1:length(g.x);
ki=sort(setdiff(k,ju));
h.x=g.x(ju);
h.y=g.y(ju);
h.z=g.z(ju);

h.e=g.e;

%for k=ki
%    j=find(g.e > k );
%    h.e(j)=h.e(j)-1;
%end

for k=1:length(ki)
    j=find(g.e > ki(k) );
    h.e(j)=h.e(j)-1;
end
