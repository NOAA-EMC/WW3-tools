function g=RemoveSandPointsWW3(g,FileOutWW3)

close all;
count=0;
IsValidBoundary=0
while ~IsValidBoundary, 
    x=g.x;y=g.y;f=g.z;e=g.e;
    bnd=detbndy(e);
    bndn=unique(bnd(:));
    clear n1 n2
    for k=1:length(bndn);
        j1=find(bndn(k)==bnd(:,1));
        j2=find(bndn(k)==bnd(:,2));
        n1(k)=length(j1);
        n2(k)=length(j2);
    end
    
    nb=n1+n2;
    unique(nb)% 2 or 4

    j4=find(nb==4)
    if isempty(j4)
        IsValidBoundary=1
    end
  
    figure;
    clf;patch(x(e'),y(e'),f(e'));shading interp;colormap('jet');colorbar;
    hold on
    plot(x(bndn(j4)),y(bndn(j4)),'ro');
    title(['iteration : ',int2str(count)])
    if ~IsValidBoundary,
        count=count+1;
         bb=bndn(j4)
         bbu=unique(bb)
         display(['number of bad boundary points= ',int2str(length(bbu))])
         display(['iteration number : ',int2str(count)])
         nn=length(x);
         %jgn=setdiff([1:nn],bbu(:)');
         
         jgnBnd=setdiff([1:nn],bbu(:)');
         g.x=x(:);g.y=y(:);g.z=f(:);g.e=e;
         hBnd=submeshFast(g,jgnBnd);%can orphan nodes 
         jgnInt=unique(hBnd.e(:));%nodes still in an element
         h=submeshFast(hBnd,jgnInt);
         [ne,three]=size(h.e)
         nn=length(h.x)
         g=h;
    end
end

if nargin >1,
    WriteWW3MeshX(g,FileOutWW3);
end