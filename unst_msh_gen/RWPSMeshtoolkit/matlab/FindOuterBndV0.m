function OpenBndNodes=FindOuterBnd(g)



e=g.tria3.index(:,1:3);
x=g.point.coord(:,1);y=g.point.coord(:,2);f=g.point.coord(:,3);

bnd=detbndy(e);
c=edges2chains(bnd);
for k=1:length(c)
    k
    xc(k)=mean(x(c(k).nodes));    
    yc(k)=mean(y(c(k).nodes));
    A(k)=areaint(x(c(k).nodes),y(c(k).nodes));
end

[m,k]=max(A);
jin=setdiff(1:length(A),k)
obj=k;
if A(obj)>sum(A(jin)),
    disp(['outer boundary is: ',int2str(obj)])
    plot(x(c(obj).nodes),y(c(obj).nodes),'r.-');
end
hold on
for k=1:length(jin)
    j=jin(k);
    plot(x(c(j).nodes),y(c(j).nodes),'k-')
end

figure;ph=patch(x(e'),y(e'),f(e'));shading interp;colormap('jet');colorbar
hold on;plot(x(c(obj).nodes),y(c(obj).nodes),'k.-');
n=c(obj).nodes;

g=loadmsh('output/RWPS.F.LLH.msh')
p=loadmsh('PSLGboundary5km_shifted.msh');
x=g.point.coord(:,1);y=g.point.coord(:,2);f=g.point.coord(:,3);
xp=p.point.coord(:,1);yp=p.point.coord(:,2);fp=p.point.coord(:,3);

%remove corners from boundary where they exist
eps=5000;%
figure;plot(xp,yp,'k.',max(xp),min(yp),'ro',min(xp),min(yp),'ro',max(xp),max(yp),'ro',min(xp),max(yp),'ro');
wgs84 = wgs84Ellipsoid("m");
DSW=distance( [min(yp),min(xp)],[yp,xp],wgs84);
DSE=distance( [min(yp),max(xp)],[yp,xp],wgs84);
DNW=distance( [max(yp),min(xp)],[yp,xp],wgs84);
DNE=distance( [max(yp),max(xp)],[yp,xp],wgs84);

jSW=find(DSW<eps);
jNW=find(DNW<eps);
jSE=find(DSE<eps);
jNE=find(DNE<eps);
np=length(xp);
jp=1:np;
jp=setdiff(jp,jSW);
jp=setdiff(jp,jNW);
jp=setdiff(jp,jSE);
jp=setdiff(jp,jNE);

xp=xp(jp);
yp=yp(jp);
figure;plot(xp,yp,'k.',max(xp),min(yp),'ro',min(xp),min(yp),'ro',max(xp),max(yp),'ro',min(xp),max(yp),'ro');

figure
plot(x(n),y(n),'ko',xp-360,yp,'r.');

%compute distance from each outer boundary point in mesh to PSLG
% defining boundary- slow cooking.
if 0,
    t0=now
    clear d
    wgs84 = wgs84Ellipsoid("m");
    for k=1:length(n)
        if mod(k,10)==0,
            k/length(n),
            etr=(length(n)-k)*(now-t0)/k;
            disp(['estimated time remaining: ',num2str(24*etr),' hours'])
        end
        D=distance( [y(n(k)),x(n(k))],[yp,xp-360],wgs84);
        d(k)=min( D );
    end
    save -v7.3 outerbounds.mat x y n d xp yp
end

%compute distance from each outer boundary point in mesh to PSLG
% defining boundary- ~100x faster.

lat2m=single(110574)
t0=now
clear d
xpp=xp-360;
for k=1:length(n)   
    lon2m=single(111320.*cos(y(n(k))*pi/180));
    if mod(k,1000)==0,
        k/length(n),
        etr=(length(n)-k)*(now-t0)/k;
        disp(['estimated time remaining: ',num2str(24*etr),' hours'])
    end
    DLON= single(mod(x(n(k))-xpp(:),360)) ;
    DLAT= single(y(n(k))-yp(:));
    ld1=min (  abs(       DLON*lon2m + i*DLAT*lat2m  ) );
    ld2=min (  abs( (360-DLON)*lon2m + i*DLAT*lat2m  ) );
    d(k)=min(ld1,ld2  );
end
save -v7.3 outerbounds.mat x y n d xp yp

epsD=3000
dx=4
joo=find(d>epsD); % points 3 km or more from coastline point
noo=n(joo);
figure;plot(xp-360,yp,'y.',x(n),y(n),'k.-',x(noo),y(noo),'r.')

xmax=max(xpp)-dx
xmin=min(xpp)+dx
ymax=max(yp)-dx
ymin=min(yp)+dx

jX=find(and( x>xmin,x<xmax ));
jY=find(and( y>ymin,y<ymax ));
jint=intersect(jX,jY);
nooint=intersect(noo,jint);
nooext=setdiff(noo,nooint);
whos noo

e=g.tria3.index(:,1:3);
clf;patch(x(e'),y(e'),f(e'));shading interp;colormap('jet');colorbar;
hold on
plot(x(nooext),y(nooext),'ro')

OpenBndNodes=nooext;