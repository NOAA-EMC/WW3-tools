function gnew=RemoveMissingIslandsEle(g,ax,p,MinBndDist,isplot);
%function gnew=RemoveMissingIslandsEle(g,ax,p,MinBndDist,isplot);
% Remove nodes of mesh, g inside closed curves defined in pslg
% p. Only nodes in which all adjacent elements are with  a closed 
% curve are removed.  The action is preformed only within the axis ax.
% This wont touch nodes within MinBndDist (m) of existing boundary
%   input:
%          g : FE mesh structure with fields
%               g.x : longitute
%               g.y : latitude
%               g.z : bathymetric depth 
%               g.e : (ne x 3) element list
%          ax : output of axis for figure, i.e. ax=axis or prespecified as ax=[xmin, xmax,ymin, ymax]
%          p: a Planar Straight Line Graph (pslg) structure with fields
%               p.x : (nn x 1) x coordinates of nodes
%               p.y : (nn x 1 )y coordinates of nodes
%               p.edges : (nedges x 2) list of edges between nodes
%          MinBndDist : distance in meters from p in which mesh nodes in g will not be removed
%          isplot=0 for no ploting and isplot=1 to plot mesh editing
%   output:
%          gnew : FE mesh structure with nodes inside curves defined by p removed
%

if nargin<4
    MinBndDist=1000;
end
if nargin <5,
    isplot=0;
end
jx=find(and(g.x>ax(1),g.x<ax(2)));
jy=find(and(g.y>ax(3),g.y<ax(4)));
jg=intersect(jx,jy);
if ~isempty(jg),
    whos
    g0=submeshFast(g,jg);
    bnd=detbndy(g0.e);
    jbnd=unique(bnd(:));

    xb=g0.x(jbnd);
    yb=g0.y(jbnd);

    jx=find(and(p.x>ax(1),p.x<ax(2)));
    jy=find(and(p.y>ax(3),p.y<ax(4)));
    jp=intersect(jx,jy);
    if ~isempty(jp),
        p0=subpslgFast(p,jp);
        x=g0.x;y=g0.y;z=g0.z;e=g0.e;
        chains=edges2chainsG(p0.edges);
        %chains=edges2chains(p0.edges);
        if isplot,
            figure(1);clf;ph=patch(x(e'),y(e'),z(e'));cm=colormap('jet');shading interp;colorbar;axis equal;
            hold on;plot(p0.x(p0.edges'),p0.y(p0.edges'),'k.-');
            set(ph,'EdgeColor','w');set(ph,'EdgeAlpha',.2);
        end

        nc=length(chains);
        jb0=[];
        for k=1:nc 
            IsClosed(k)=0;
            n=chains(k).nodes;
            if length(n)>2,
                if n(1)==n(end),
                    %jins=find(inside(x,y,p0.x(n),p0.y(n)));
                    %jins=find(inside(x,y,p0.x(n),p0.y(n))==1);
                    [inP,onP]=insidepoly(x,y,p0.x(n),p0.y(n));
                    jins=find(or(inP,onP));

                    jb0=[jb0;jins(:)];
                    IsClosed(k)=1;
                    hold on;plot(p0.x(n),p0.y(n),'gx-')
                end
            end
        end

        lat2m=110574.;
        lon2m=111320.;

        jb0X=[];
        for k=1:length(jb0);
            n=jb0(k);

            d=min(abs( lon2m*( x(n)-xb )*cos(pi*y(n)/180) + ...
                    i*lat2m*( y(n)-yb ) ));
            if d>MinBndDist,jb0X=[jb0X,n];end
        end
        jb0=jb0X;

        if isplot,
            hold on
            plot(x(jb0),y(jb0),'ro');
        end

        [g1,kr0]=submeshFastEle(g0,jb0);

        if isplot,
            figure(2);
            x=g1.x;y=g1.y;z=g1.z;e=g1.e;
            clf;ph=patch(x(e'),y(e'),z(e'));cm=colormap('jet');shading interp;colorbar;axis equal;
            hold on;plot(p0.x(p0.edges'),p0.y(p0.edges'),'k.-');
            set(ph,'EdgeColor','w');set(ph,'EdgeAlpha',.2);
        end
        %jb=jg(jb0);
        kr1=intersect(kr0,jb0);% kr0 can have boundary nodes not intended for removal!!
        if isplot,
            figure(2);
            XE=g0.x(kr1);
            YE=g0.y(kr1);
            plot(XE,YE,'md')
            figure(1);hold on;plot(XE,YE,'md')
        end
        jb=jg(kr1);
        nn=length(g.x);
        gnew=submeshFast(g,setdiff(1:nn,jb)); %only removes nodes
    else
        gnew=g;
    end
else
    gnew=g;
end

