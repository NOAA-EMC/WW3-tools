function [xss,yss]=SmoothSubSampleCoastlineFast(x,y,DI,lambda);
%function [xs,ys]=SmoothSubSampleCoastline(x,y,dsmooth,lambda);
%Sub sample coastline at dsmooth (m) distance. Boxcar smooth coastline to 
% dsmooth*lambda distance   
itz=0;
[xi,yi]=InterpCoastline1(x,y,DI,itz);
gamma=2*lambda+1;
ic=lambda+1;
W=ones(1,gamma)/gamma;
xs=conv(xi,W,'same');
ys=conv(yi,W,'same');
xss=xs(ic:lambda:end-ic);
yss=ys(ic:lambda:end-ic);
if and(x(1)==x(end),y(1)==y(end))
    if length(xss)>0,
        xss=[xss,xss(1)];
        yss=[yss,yss(1)];
    end
end
