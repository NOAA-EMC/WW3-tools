function PlotWW3Mesh(flnm,ax,caxZ,caxLS);
% function PlotWW3Mesh(flnm,ax,caxZ,caxLS)
%
% Plot figure showing colorfields for bathymetry and 
% lengthscale for an unstructured mesh in WW3 
% file format. Exports figures in .jpg files in a 
% "-nodesktop" environment.
%
%   inputs:
%       flnm : name of file defining an unstructured 
%              mesh in WW3 .msh format
%       ax : [West,East,South,North] window to make a zoom 
%            plot. this is in addition to a figure showing
%            the full domain.
%       caxZ(optional) : [minVal,maxVal] to limit bathymetry 
%                        colorfield
%       caxLS(optional) : [minVal,maxVal] to limit length scale 
%                        colorfield
%
% example call:
% >> PlotWW3Mesh('RWPS.WW3a1.msh',[-123.0885 -122.6996   48.3679   48.7568],[0,3000],[0,8]);
%


res=300;% resolution to print graphics with 

set(groot, 'defaultFigureVisible', 'off'); 
figure; 

g=loadmshWW3(flnm)
x=g.x;y=g.y;z=g.z;e=g.e;
x=LonCon(x)

clf;
ph=patch(x(e'),y(e'),z(e'));
shading interp;
colorbar('v')
cm=colormap('jet');colormap(flip(cm));

if nargin>2,caxis(caxZ);end

daspect([1,cos(mean(y)*pi/180),1]);
title('mesh bathymetry(m)')
exportgraphics(gcf, [flnm(1:end-4),'.bathy.jpg'], 'Resolution', res);

if nargin>1
    axis(ax)
    exportgraphics(gcf,[flnm(1:end-4),'.ax.bathy.jpg'],'Resolution', res);
    set(ph,'EdgeColor','k');set(ph,'EdgeAlpha',.2);
    exportgraphics(gcf,[flnm(1:end-4),'.edges.ax.bathy.jpg'],'Resolution', res);
end

LS=ComputeLengthScale_wgs84_MEL(x,y,e);
LSn=Ele2Nodes(x,y,e,LS);
clf;
ph=patch(x(e'),y(e'),LSn(e'));
shading interp;
colorbar('v')
cm=colormap('jet');colormap(flip(cm));

if nargin>3,caxis(caxLS);end

daspect([1,cos(mean(y)*pi/180),1]);
title('mesh length sclae (km)')
exportgraphics(gcf,[flnm(1:end-4),'.lengthscale.jpg'],'Resolution', res);

if nargin>1
    axis(ax)
    exportgraphics(gcf,[flnm(1:end-4),'.ax.lengthscale.jpg'],'Resolution', res);
    set(ph,'EdgeColor','k');set(ph,'EdgeAlpha',.2);
    exportgraphics(gcf,[flnm(1:end-4),'.edges.ax.lengthscale.jpg'],'Resolution', res);
end



close(gcf);

