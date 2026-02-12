function  Fnode = Ele2Nodes(lon,lat,e,Fele)
%function  Fnode = Ele2Nodes(lon,lat,e,Fele)
% Approximate element averages for field supported on elements 
% at nodes.  Uses average weighted by elment areas
% inputs:
%       lon : (nn x 1) longitude coordinates of nodes
%       lat : (nn x 1) latitude coordinates of nodes
%       e   : (ne x 3) element matrix
%       Fele   : (ne x 1) field average on elements
% output:
%       Fnode : (nn x 1) approximation of Fele at nodes 
%

A=EleArea(lon,lat,e);
nn=length(lon);
[ne,three]=size(e);
Fnode=zeros(nn,1);
Anode=zeros(nn,1);
for k=1:ne
    if mod(k,100000)==0,k/ne,end
    j=e(k,:);
    Fnode(j)=Fnode(j)+A(k)*Fele(k);
    Anode(j)=Anode(j)+A(k);
end
Fnode=Fnode./Anode;
