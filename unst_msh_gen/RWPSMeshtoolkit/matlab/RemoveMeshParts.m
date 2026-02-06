function gnew=RemoveMeshParts(g,ax,S,ca);

%function gnew=HandEditMesh(g,ax,p);
% Graphical interface to remove parts of mesh, g inside axis, where ax=[xmin, xmax,ymin, ymax].
% The frame work is intended to expand to moving nodes, removing elements, etc
%
%   inputs:
%           g : FE mesh structure with fields
%               g.x : longitute
%               g.y : latitude
%               g.z : bathymetric depth 
%               g.e : (ne x 3) element list
%           ax : output of axis for figure, i.e. ax=axis or prespecified as ax=[xmin, xmax,ymin, ymax]
%           S : structure output from S=shaperead(CoastLineFile) etc.
%           ca : [1,2] caxis for colorfield based on g.z ca=[min(g.z),max(g.z)]
%

jx=find(and(g.x>ax(1),g.x<ax(2)));
jy=find(and(g.y>ax(3),g.y<ax(4)));
jg=intersect(jx,jy);
g0=submeshFast(g,jg);
bnd=detbndy(g0.e);
jbnd=unique(bnd(:));

xb=g0.x(jbnd);
yb=g0.y(jbnd);

x=g0.x;y=g0.y;z=g0.z;e=g0.e;

figure;
clf;ph=patch(x(e'),y(e'),z(e'));cm=colormap('jet');shading interp;colorbar;axis equal;caxis(ca);axis(ax)
hold on;
for k=1:length(S),plot(S(k).X,S(k).Y,'k');end
 set(ph,'EdgeColor','w');set(ph,'EdgeAlpha',.2);
 axis(ax)


b=1;
rmi=[]
while b==1
    disp('(left button) Drag and drop a rectangle over region to delete- right click to quit');
    title('(left button) Drag and drop a rectangle over region to delete- right click to quit');
    [x0,y0,x1,y1,b]=dragdrop;
    if b==1,
        s=sum([g0.x(:)'>x0;g0.x(:)'<x1;g0.y(:)'>y0;g0.y(:)'<y1]);
        bind=find(s==4)%chuckem nodes
        hold on;
        plot(g0.x(bind),g0.y(bind),'rx');
        rmi=[rmi,bind];
    end
end

title('Now click on individual nodes to remove, right click when done');
bb=1

while bb==1,
    [xx,yy,bb]=ginput(1)
    [mm,jj]=min(abs(g0.x+i*g0.y-xx-i*yy));
    if bb==1
        rmi=[rmi,jj];
        plot(g0.x(jj),g0.y(jj),'kx')
    end
end

nn=length(g0.x);
jg0=setdiff(1:nn,rmi);
g1=submeshFast(g0,jg0);

figure;
x=g1.x;y=g1.y;z=g1.z;e=g1.e;
clf;ph=patch(x(e'),y(e'),z(e'));cm=colormap('jet');shading interp;colorbar;axis equal;caxis(ca)
hold on;
for k=1:length(S),plot(S(k).X,S(k).Y,'k');end
anodes=unique(e(:));
nn=length(x);
deadnodes=setdiff(nn,anodes);
if length(deadnodes)>0,
    hold on;plot(x(deadnodes),y(deadnodes),'ko',x(deadnodes),y(deadnodes),'kx');
    title('Yo! nodes now exist with no elements! Fix with remove_dead_nodes');
end

set(ph,'EdgeColor','w');set(ph,'EdgeAlpha',.2);
axis(ax)
whos jg rmi
max(rmi)
jb=jg(rmi);
nn=length(g.x);
gnew=submeshFast(g,setdiff(1:nn,jb));


