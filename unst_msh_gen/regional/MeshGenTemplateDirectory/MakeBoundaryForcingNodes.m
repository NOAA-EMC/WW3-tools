function MakeBoundaryForcingNodes(dx,CheckGFS)
if nargin<2,CheckGFS=0;end

Blon=[129.91 10.71];
Blon(1)=Blon(1)-360;
Blat=[-30.42 79.99];
g=loadmshWW3('RWPS.OSMxGSHHS.GMU.txt.WW3.msh');

jb=unique(g.bnd);

xb=g.x(jb);
yb=g.y(jb);

%x=Blon(1)-dx:dx:Blon(2)+dx;
%y=Blat(1)-dx:dx:Blat(2)+dx;
x=Blon(1):dx:Blon(2);
y=Blat(1):dx:Blat(2);
if x(end)~=Blon(2),
    x=[x,Blon(2)];
end
if y(end)~=Blat(2),
    y=[y,Blat(2)];
end


Xs=x;
Ys=zeros(size(Xs))+min(Blat);
Xn=x(end:-1:1);
Yn=zeros(size(Xs))+max(Blat);

Ye=y
Xe=zeros(size(Ye))+max(Blon);
Yw=y(end:-1:1);
Xw=zeros(size(Yw))+min(Blon);

X=[Xs,Xe,Xn,Xw];
Y=[Ys,Ye,Yn,Yw];
whos X Y
clear m i;
for k=1:length(X);
    m(k)=min(abs(X(k)+i*Y(k)-xb-i*yb));
end

j=find(m<2*dx);
Xp=X(j);
Yp=Y(j);
Xpp=Xp;
jW=find(Xpp<0);
Xpp(jW)=Xpp(jW)+360;

filename = ['BoundaryDX',int2str(round(1/dx)),'thDeg.RWPS.txt']
fileID = fopen(filename, 'w');
fprintf(fileID, '%6f %6f\n',[Xpp(:),Yp(:)]');
fclose(fileID);

figure;
x=g.x;y=g.y;z=g.z;e=g.e;
clf;ph=patch(x(e'),y(e'),z(e'));shading interp;colormap('jet');colorbar; caxis([0,5000]);axis equal
hold on
plot(Xp,Yp,'k.-');
plot(xb,yb,'r.')

if CheckGFS,
    Xpp=Xpp(:);
    Yp=Yp(:);
    xb=Xpp;
    yb=Yp;
     g=loadmshWW3('/mnt/sda/keston/meshes/uglo_15km.msh');
    x=g.x;y=g.y;z=g.z;e=g.e;
    
    clf;patch_global(x,y,z,e');shading interp;
    cm=colormap('jet');colormap(flip(cm));colorbar; caxis([0,5000]);axis equal
    
    zb=scattered_interp2_noext(g.x,g.y,g.z,g.e,xb,yb)
    j=find(isnan(zb ))
    
    nb=length(xb);
    k=setdiff(1:nb,j)
    hold on;
    phg=plot(xb(k),yb(k),'w.');
    phb=plot(xb(j),yb(j),'ko',xb(j),yb(j),'kx');
    set(phb(1),'MarkerSize',13);set(phb(2),'MarkerSize',13);
    title([int2str(round(1/dx)),'th degree spacing. White dot RWPS bound in GFS mesh, Black-Circle points out of GFS mesh ']);
    kprint(['BoundaryDX',int2str(round(1/dx)),'thDeg.RWPS.jpg'])
    InGFS=k;
    OutGFS=j
     
    filename = ['BoundaryDX',int2str(round(1/dx)),'thDeg.RWPS.InGFS.txt']
    fileID = fopen(filename, 'w');
    fprintf(fileID, '%6f %6f\n',[Xpp(k),Yp(k)]');
    fclose(fileID);
end
