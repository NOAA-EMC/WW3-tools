function j=FindPointsAx(ax,x,y);
%function j=FindPointsAx(ax,x,y);
% returns indexs j for points (x,y) inside matlab axis object ax
% useage example: 
%
% >>ax=axis % get current figure axis
% >>j=FindPointsAx(ax,g.x,g.y);
% >>g_local=submeshFast(g,j)
% 
% to return mesh structure contained within current figure axis

jx=find(and(x>ax(1),x<ax(2)));
jy=find(and(y>ax(3),y<ax(4)));
j=intersect(jx,jy);
