function j=FindPointsAx(ax,x,y);
jx=find(and(x>ax(1),x<ax(2)));
jy=find(and(y>ax(3),y<ax(4)));
j=intersect(jx,jy);