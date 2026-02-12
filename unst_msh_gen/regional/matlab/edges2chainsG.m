function chains=edges2chainsG(edges);
%function chains=edges2chainsG(edges);
% This performs similair function to edges2chains, how ever it uses matlabs graph theory library
% and only returns closed cycles rather than cycles and open ended chains. 
% Somewhat faster than origonal edges2chains.m  (~x 2 speedup)
%

G=graph(edges(:,1),edges(:,2));
cycles = allcycles(G);
cc=0;
for k=1:length(cycles),
    n=cycles{k};
    if length(n)>2,
        n=[n(:);n(1)];
        cc=cc+1;
        chains(cc).nodes=n;
    end
end
