function h=BoundingBox(x,y,c);
%function h=BoundingBox(x,y,c);
% draw a bounding rectanble around points x, y in color c. 
% Note, x and y can be different lengths here.

if nargin<3,c='k';end
x0=min(min(x));
y0=min(min(y));
x1=max(max(x));
y1=max(max(y));

hold on;
h=plot([x0,x1,x1,x1,x0,x0],...
     [y0,y0,y1,y1,y1,y0],c);

