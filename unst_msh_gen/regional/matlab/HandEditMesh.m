function gnew=HandEditMesh(g,ax,p);
%function gnew=HandEditMesh(g,ax,p);
% Hand remove parts of mesh, g, inside closed curves defined in pslg p.
% The action is preformed only within the axis ax=[min x, max x, min y , max y].
% The frame work is intended to expand to moving nodes, removing elements, etc


jx=find(and(g.x>ax(1),g.x<ax(2)));
jy=find(and(g.y>ax(3),g.y<ax(4)));
jg=intersect(jx,jy);
g0=submeshFast(g,jg);
bnd=detbndy(g0.e);
jbnd=unique(bnd(:));

xb=g0.x(jbnd);
yb=g0.y(jbnd);

jx=find(and(p.x>ax(1),p.x<ax(2)));
jy=find(and(p.y>ax(3),p.y<ax(4)));
jp=intersect(jx,jy);
p0=subpslgFast(p,jp);

x=g0.x;y=g0.y;z=g0.z;e=g0.e;
chains=edges2chains(p0.edges);

clf;ph=patch(x(e'),y(e'),z(e'));cm=colormap('jet');shading interp;colorbar;axis equal;
hold on;plot(p0.x(p0.edges'),p0.y(p0.edges'),'k.-')
set(ph,'EdgeColor','w');set(ph,'EdgeAlpha',.2);

title('left click on mesh nodes to remove, right click when done')
bind=[];
b=1
while b==1,
	[xx,yy,b]=ginput(1);
   if b==1,
      [m,ind]=min(abs(xx+i*yy-x-i*y));
   	plot(x(ind),y(ind),'rx');   
   	bind=[bind,ind];   
   end
end

jb0=bind;
nn=length(x);
jg0=setdiff(1:nn,jb0);
g1=submeshFast(g0,jg0);

figure;
x=g1.x;y=g1.y;z=g1.z;e=g1.e;
clf;ph=patch(x(e'),y(e'),z(e'));cm=colormap('jet');shading interp;colorbar;axis equal;
hold on;plot(p0.x(p0.edges'),p0.y(p0.edges'),'k.-')
set(ph,'EdgeColor','w');set(ph,'EdgeAlpha',.2);

jb=jg(jb0);
nn=length(g.x);
gnew=submeshFast(g,setdiff(1:nn,jb));


