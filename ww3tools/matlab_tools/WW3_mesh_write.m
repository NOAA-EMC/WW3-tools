function [filename]=WW3_mesh_write(tri,x,y,h,OB_ID,filename,plott)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function writes WW3 grid including nodes             %
% (longitude,latitude,depth), open bounday nodes            %
% and element connections (triangles)                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Ali Abdolali August 2018 ali.abdolali@noaa.gov       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tri: triangle (nelem,3);
% x: longitude (nnode,1);
% y: laritude (nnode,1);
% h: depth (nnode,1);
% OB_ID: Open boundary nodes ID
% filename: mesh name
% Plot it, if plott=1

node(:,1)=(1:length(x));
node(:,2)=x;
node(:,3)=y;
node(:,4)=h;

fileID = fopen(filename,'w');
fprintf(fileID,'%s\n', '$MeshFormat');
fprintf(fileID,'%s\n', '2 0 8');
fprintf(fileID,'%s\n', '$EndMeshFormat');
fprintf(fileID,'%s\n', '$Nodes');
fprintf(fileID,'%d\n', length(node(:,1)));

for i=1:length(node(:,1))
    fprintf(fileID,['%d %s %5.5f %s %5.5f %s %5.5f\n'], node(i,1),'', node(i,2),'',node(i,3),'',node(i,4));
end
fprintf(fileID,'%s\n', '$EndNodes');
fprintf(fileID,'%s\n', '$Elements');
fprintf(fileID,'%d\n', length(tri(:,1))+length(OB_ID));
m=0;
for i=1:length(OB_ID)
    m=m+1;
    fprintf(fileID,['%d %s %d %s %d %s %d %s %d %s %d\n'], m,'',15,'',2,'',0,'',0,'',OB_ID(i));
end

for i=1:length(tri(:,1))
    m=m+1;
    fprintf(fileID,['%d %s %d %s %d %s %d %s %d %s %d %s %d %s %d %s %d\n'], m,'',2,'',3,'',0,'',i,'',0,'',tri(i,1),'',tri(i,2),'',tri(i,3));
end
fprintf(fileID,'%s', '$EndElements');
fclose(fileID);

if plott==1
width=880;  % Width of figure for movie [pixels]
height=700;  % Height of figure of movie [pixels]
left=700;     % Left margin between figure and screen edge [pixels]
bottom=200;  % Bottom margin between figure and screen edge [pixels]
 


    figure
set(gcf,'Position', [left bottom width height])



trisurf(tri,node(:,2),node(:,3),node(:,4));
shading interp

view(2);
axis equal
hold on
p1=scatter(node(OB_ID,2),node(OB_ID,3),'xk')
hCbar=colorbar
colormap(jet)
legend([p1],'Open Boundary Nodes')
xlim([min(node(:,2))-(max(node(:,2))-min(node(:,2)))/10 max(node(:,2))+(max(node(:,2))-min(node(:,2)))/10])
ylim([min(node(:,3))-(max(node(:,3))-min(node(:,3)))/10 max(node(:,3))+(max(node(:,3))-min(node(:,3)))/10])

if min(node(:,4))~=max(node(:,4))
caxis([min(node(:,4)) max(node(:,4))])
end

hold on
xlabel('Longitude ^{\circ} ');
ylabel('Latitude ^{\circ} ');
axis equal
box on
grid off
axis on
title(['WW3 Grid - depth [m]'],'fontsize',15)

print(gcf,'-dpng',[filename,'.png'],'-r900');
end
