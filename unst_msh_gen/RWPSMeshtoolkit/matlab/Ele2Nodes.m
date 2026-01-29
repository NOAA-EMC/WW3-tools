function  Fnode = Ele2Nodes(lon,lat,e,Fele)
%crude compute element Area

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
