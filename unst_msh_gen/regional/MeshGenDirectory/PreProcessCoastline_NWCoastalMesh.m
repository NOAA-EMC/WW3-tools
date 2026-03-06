
SetPath

lonWest=-130;lonEast=-121;
latSouth=45;latNorth=51;
CoastLineFile = 'GlobalCoastlineOSM.shp'

MakeCoastalBoundariesOSM
% Use SmoothCoastalBoundaries (rather than script MakeCoastalBoundariesOSM) to set
% the lengthscale of the coastline.  The call below smooths and resamples the 
% coastline to 250m rather than 500m.
% SmoothCoastalBoundaries(GlobalCoastlineFileOSM,0.25,CoastLineFile) 

BuildBoundaryPSLGfunction(CoastLineFile,lonWest,lonEast,latSouth,latNorth,PSLGfile)
%[xpi,ypi]=ExtraPointsOfIntrest;

TargetShape=TargetCoastlineFile;% US coastline as target for distance
DX=.0125;
MakeDistanceToCoastData(DX,TargetCoastlineFile,[],[])
