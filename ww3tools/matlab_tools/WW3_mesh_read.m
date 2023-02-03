function [tri,x,y,z] = WW3_mesh_write(file,plott)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function reads WW3 grid including nodes              %
% (longitude,latitude,depth) and element connections        %
%  (triangles)                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Ali Abdolali August 2018 ali.abdolali@noaa.gov       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if plot==1, then plot the mesh

% no of nodes is mentioned in 5th row and first column

N_n      = dlmread(file,'',[5-1 1-1 5-1 1-1]);
N_e      = dlmread(file,'',[7+N_n 0 7+N_n 0]);

node_id     = dlmread(file,'',[5 0 4+N_n 0]);
element_id     = dlmread(file,'',[8+N_n 0 8+N_e+N_n-1 0]);

nodes       = dlmread(file,'',[5 1 4+N_n 3]);
TRI    = dlmread(file,'',[8+N_n 0 8+N_e+N_n-1 8]);

%------- 2D Geometry

two_d_nodes = nodes(:,1:2);
elem_type   = TRI(:,3);

%--- find the starting indices of 2D elements
two_ind = 1;
for i = 1:N_e
    if(elem_type(i) ~= 2)
        two_ind = two_ind+1;
    end
end
%----------------------------------------------

%two_d_elements(1:N_e-two_ind,1:3) = 0;
   two_d_elements(:,1:3) = TRI(N_e-two_ind+2:N_e,7:9);
   
   tri=two_d_elements;
   x=nodes(:,1);
   y=nodes(:,2);
   z=nodes(:,3);
   
   
   %%
   
   if plott==1;
width=880;  % Width of figure for movie [pixels]
height=700;  % Height of figure of movie [pixels]
left=700;     % Left margin between figure and screen edge [pixels]
bottom=200;  % Bottom margin between figure and screen edge [pixels]

figure
set(gcf,'Position', [left bottom width height])


cmap = colormap;

trisurf(tri,x,y,z);
shading interp
view(2)
hcb=colorbar
caxis([min(z) max(z)])
title(hcb,'[m]','fontsize',14);
colormap(jet)
hold on

hold on
 
sz=40;

 
xlabel('Longitude ^o ');
ylabel('Latitude ^o ');

title(['depth (m)'],'fontsize',15)



ylim([min(y)-0.1*((max(y)-min(y))) max(y)+0.1*((max(y)-min(y)))]);
xlim([min(x)-0.1*((max(x)-min(x))) max(x)+0.1*((max(x)-min(x)))]);

grid off
box on


 print('-dpng',[file,'.png']);
   end
 