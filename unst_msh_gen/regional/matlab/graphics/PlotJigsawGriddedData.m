function PlotJigsawGriddedData(flnm,ax,cax);
% function PlotJigsawGriddedData(flnm,ax,cax)
%
% Plot figure showing colorfields for gridded data 
% in jigsaw file format. Exports figures in .jpg 
% files in a "-nodesktop" environment.
%
%   inputs:
%       flnm : name of file defining an gridd data
%              in jigsaw format
%       ax : [West,East,South,North] window to make a zoom 
%            plot. this is in addition to a figure showing
%            the full domain.
%       cax(optional) : [minVal,maxVal] to limit colorfield
%
% example call:
% >> PlotJigsawGriddedData('DFun.NWcoastal.PSLG.msh',[-123.0885 -122.6996   48.3679   48.7568],[0,3*10^5]);
% >> PlotJigsawGriddedData('Topo.DFun.NWcoastal.PSLG.msh',[-123.0885 -122.6996   48.3679   48.7568],[0,300]);
%

res=300;% resolution to print graphics with 

set(groot, 'defaultFigureVisible', 'off'); 
figure; 

g=loadmsh(flnm)
x=g.point.coord{:,1};
y=g.point.coord{:,2};
z=g.value;

clf;
ph=pcolor(x,y,z);
shading interp;
colorbar('v')
cm=colormap('jet');%colormap(flip(cm));

if nargin>2,caxis(cax);end

daspect([1,cos(mean(y)*pi/180),1]);
title(['field from ',flnm])
exportgraphics(gcf, [flnm(1:end-4),'.field.jpg'], 'Resolution', res);

if nargin>1
    axis(ax)
    daspect([1,cos(mean(y)*pi/180),1]);
    exportgraphics(gcf,[flnm(1:end-4),'.ax.field.jpg'],'Resolution', res);
    set(ph,'EdgeColor','k');set(ph,'EdgeAlpha',.2);
    exportgraphics(gcf,[flnm(1:end-4),'.ax.edges.field.jpg'],'Resolution', res);
end

close(gcf);
