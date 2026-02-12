function [xss,yss]=SmoothSubSampleCoastlineFast(x,y,DI,lambda);
% function [xss,yss]=SmoothSubSampleCoastline(x,y,dsmooth,lambda);
% Smooth and subsample coastline at dsmooth (m) distance. The coastline is first interpolated
% to uniform spacing and then a boxcar smoother of width DI*(2*lambda+1)is applied.
%   Inputs:
%       x,y :   longitude, latitude points along a coastline. 
%               If the coastline defines an island 
%               then x(1)==x(end) and y(1)==y(end),
%       DI: distance in meters to interpolate the coastline to
%       lambda : integer coastline is smoothed to DI*(2*lambda+1) and resampled lambda*DI spacing
%   outputs : 
%       xss,yss : longitude, latitude points the smoothed and subsampled coastline
%
itz=0;
[xi,yi]=InterpCoastline(x,y,DI,itz);
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
