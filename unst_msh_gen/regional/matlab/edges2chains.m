function chains=edges2chains(edges);
%edges2chains::decompose list of edges into all maximal chains and cycles(non edge repeating paths)
%%take a Nx2 adjacency list and decompose it into all maximum cylcles and chains with provided no 
%vertex has order greater than 2.
%function chains=edges2chains(edges);
%
%Inputs:
%       edges:(Nedgesx2) integer adjacency list defining edges in graph.
%
%Outputs:
%       chains:list of structures defining maximal paths.
%           chains(k).nodes:n_k length integer list of nodes in chain k.
%
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%Keston Smith & Ata Bilgili May 2001
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


if ~isempty(edges)
    edges=[min(edges')',max(edges')'];
    [i,j]=sort(edges(:,1));

    edges=[edges(j,1),edges(j,2)];

    [N,NN]=size(edges);
    chains=[];
    n=0;
    while ~isempty(edges)
        chain=edges(1,1);
        j=edges(1,2);
        edges=edges(2:end,:);
        %while ~isempty(j),
        while length(j)==1,
            chain=[chain,j(:)'];
            %whos edges j
            k1=find(edges(:,1)==j);
            k2=find(edges(:,2)==j);
            if ~isempty(k1),
                j=edges(k1,2);
                [N,M]=size(edges);
                edges=edges(setdiff(1:N,k1),:);
           elseif ~isempty(k2),
                j=edges(k2,1);
                [N,M]=size(edges);
                edges=edges(setdiff(1:N,k2),:);
            else 
               j=[];        
            end
        end
        j=chain(1);
        k1=find(edges(:,1)==j);
        k2=find(edges(:,2)==j);
        if ~isempty(k1),
            j=edges(k1,2);
            [N,M]=size(edges);
            edges=edges(setdiff(1:N,k1),:);
        elseif ~isempty(k2),
            j=edges(k2,1);
            [N,M]=size(edges);
            edges=edges(setdiff(1:N,k2),:);
        else 
            j=[];        
        end
%        while ~isempty(j),
        while length(j)==1,
   
            chain=[j,chain];
            k1=find(edges(:,1)==j);
            k2=find(edges(:,2)==j);
            if ~isempty(k1),
                j=edges(k1,2);
                [N,M]=size(edges);
                edges=edges(setdiff(1:N,k1),:);
            elseif ~isempty(k2),
                 j=edges(k2,1);
                 [N,M]=size(edges);
                 edges=edges(setdiff(1:N,k2),:);
            else 
               j=[];        
            end
        end
    n=n+1;
    chains(n).nodes=chain;
    end
else
    chains(1).nodes=[];
end
