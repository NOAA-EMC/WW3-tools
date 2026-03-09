function PlotJigsawPSLG(flnm,ax)
% function PlotJigsawPSLG(flnm,ax)
%
% Plot figure showing Piecewise Straight Line Graph(PSLG)
% and print to a .jpg file in a "-nodesktop" environment
%
%   inputs:
%       flnm- name of file defining a PSLG in jigsaw format
%       ax=[West,East,South,North] window to make a zoom 
%           plot. this is in addition to a figure showing
%           the full domain.
%
% example:
% >>PlotJigsawPSLG('NWcoastal.PSLG.msh',[-123.0885 -122.6996   48.3679   48.7568]);
%

res=600;% resolution to print graphics with 

set(groot, 'defaultFigureVisible', 'off'); 
figure; 

g=loadmsh(flnm)
x=g.point.coord(:,1);
y=g.point.coord(:,2);
edg=g.edge2.index(:,1:2)';

clf;
ph=plot(x,y,'k.',x(edg),y(edg),'r');
daspect([1,cos(mean(y)*pi/180),1]);

title(['PSLG from ',flnm])
exportgraphics(gcf, [flnm(1:end-4),'.jpg'], 'Resolution', res);

if nargin>1
    axis(ax)
    daspect([1,cos(mean(y)*pi/180),1]);
    exportgraphics(gcf,[flnm(1:end-4),'.ax.jpg'],'Resolution', res);
end

close(gcf);
