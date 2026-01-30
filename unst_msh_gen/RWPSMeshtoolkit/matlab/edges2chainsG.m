function chains=edges2chainsG(edges);

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
