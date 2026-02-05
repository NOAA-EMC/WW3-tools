function pslg=joinpslg(pslg1,pslg2,mergeDist);
% function pslg=joinpslg(pslg1,pslg2,mergeDist);
% mergepslg::join two pslgs together.
%   inputs: 
%       pslg1,pslg2 : pslg structures with form
%       pslg.x : [n] point list of longitude
%       pslg.y : [n] point list of latitude
%       pslg.edges : [nedge x 2] list of edges in pslg (point adjacency)
%%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%Keston Smith 2022
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
pslg1.x=pslg1.x(:)';
pslg1.y=pslg1.y(:)';
pslg2.x=pslg2.x(:)';
pslg2.y=pslg2.y(:)';


np1=max(max(pslg1.edges));
%pslg1,pslg2,np1
if ~isempty(np1),
    pslg.edges=[pslg1.edges;pslg2.edges+np1];
else
    pslg.edges=[pslg1.edges;pslg2.edges];

end
pslg.x=[pslg1.x,pslg2.x];
pslg.y=[pslg1.y,pslg2.y];
holes.x=[];holes.y=[];
if isfield(pslg1,'holes')
    holes.x=[holes.x(:);pslg1.holes.x(:)];
    holes.y=[holes.y(:);pslg1.holes.y(:)];
end
if isfield(pslg2,'holes')
    holes.x=[holes.x(:);pslg2.holes.x(:)];
    holes.y=[holes.y(:);pslg2.holes.y(:)];
end
pslg.holes=holes;

zones.x=[];zones.y=[];
if isfield(pslg1,'zones')
    zones.x=[zones.x(:);pslg1.zones.x(:)];
    zones.y=[zones.y(:);pslg1.zones.y(:)];
end
if isfield(pslg2,'zones')
    zones.x=[zones.x(:);pslg2.zones.x(:)];
    zones.y=[zones.y(:);pslg2.zones.y(:)];
end
pslg.zones=zones;

if isfield(pslg1,'chains')
    pslg.chains=pslg1.chains;
end

if isfield(pslg2,'chains')
   nc=length(pslg.chains);
   nc2=length(pslg2.chains);
   np=length(pslg.x);
   for k=1:nc2
     pslg.chains(nc+k).nodes=n+pslg2.chains(k).nodes;
   end
end

%z=round((pslg.x+i*pslg.y)/MergeDist);
if nargin>2
    z=pslg.x+i*pslg.y;
    n=length(pslg.x);
    for k=1:n-1,
        d=abs(z(k)-z(k+1:n));
        j=find(d<mergeDist);
        pslg.edges(find(pslg.edges==j))=k;
        pslg=subpslg(pslg,j);
    end
end

