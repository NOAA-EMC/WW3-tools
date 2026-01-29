

mesh= loadmsh('mesh_mixed_coastline_planar_ocean.msh');
Elm=compute_area_Ell(mesh)
h=PatchGlobalE(mesh,Elm/1000);

 caxis([0,20]);
 colorbar
 axis([-180,180,-80,80]);
 title('edge length (km)');
 colormap(jet(20))
 kprint('EdgeLengthGlobal.jpg')
hold on

if(0)
    S=shaperead('GlobalCoastline1kmUSto15km.shp');
    N=length(S);
    xt=[];
    yt=[];
    for k=1:N
        x=S(k).X;
        y=S(k).Y;
        xt=[xt;x(:)];
        yt=[yt;y(:)];
        if mod(k,1000)==0,k/N,end
    end
    clear S
    save -v7.3 CoastVec1kmUSto15km.mat xt yt
else
    load CoastVec1kmUSto15km.mat
end
nms={'CapeCod', 'MidAtlantic', 'OlympicP','Hawaii','Spain'}
flon = [-70.6,-75,-123.5,-155.6 -2.2] ;  
flat = [41.6,38.7,47.6,  20, 41.4];
title('edge length (km)')
dx=[.25,.5,.5,.5,5]
jnan=find(isnan(xt+yt));
lx=[10,20,20,20,50]
cx=[1.5,1.5,1,1.5,30]
for k=1:5,
    dx0=dx(k)
    ax=[flon(k)-dx0,flon(k)+dx0,flat(k)-dx0,flat(k)+dx0];
     axis(ax);
     set(h,'EdgeAlpha',.1); set(h,'EdgeColor','k');
  legx=ax(1)+dx0/10;
    legy=ax(4)-dx0/10;
    [ux,uy,zn]=ll2utm(legy,legx);
    [legy5,legx5]=utm2ll(ux+lx(k)*1000,uy,zn);
    hl=plot([legx,legx5],[legy,legy5],'ro-')
    ht=text([legx+legx5]/2,legy-.05,[int2str(lx(k)),' km'])
    jx=find(and( xt<ax(2),xt>ax(1)));
    jy=find(and( yt<ax(4),yt>ax(3)));
    j=intersect(jx,jy);
    j=union(j,jnan);
    
     hold on
    plot(xt(j),yt(j),'k-');
    title([nms{k},' edge length (km)']);
    caxis([0,cx(k)]);
    colormap(jet(20));
    daspect([1,cos(flat(k)*pi/180),1])
    kprint([nms{k},'EdgeLengthGlobal.jpg']);
end

