function OpenBndNodes=FixOuterBndSandPoints(outdir)

g=loadmsh([outdir,'/RWPS.F.LLH.msh']);

e=g.tria3.index(:,1:3);
x=g.point.coord(:,1);y=g.point.coord(:,2);f=g.point.coord(:,3);

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
clf;patch(x(e'),y(e'),f(e'));shading interp;colormap('jet');colorbar;
hold on
plot(x(bndn(j4)),y(bndn(j4)),'ro');

 bb=bndn(j4)
 bbu=unique(bb)
 nn=length(x);
 gn=setdiff([1:nn]',bbu);
 h=submesh(g,gn);
save -v7.3 hmesh.mat h

[ne,three]=size(h.e)
nn=length(h.x)
g0=g;
g0=rmfield(g0,"tria3")
g0=rmfield(g0,"point")
ep=[h.e,zeros(ne,1)];
g0.tria3.index=ep;
g0.point.coord=[h.x(:),h.y(:),h.z(:),zeros(nn,1)];
g0.value=-h.z;

savemsh([outdir,'/RWPS.F.LLH.NSP.msh'],g0);
WriteWW3Mesh(g0,[outdir,'/RWPS.F.LLH.NSP.ww3'])
