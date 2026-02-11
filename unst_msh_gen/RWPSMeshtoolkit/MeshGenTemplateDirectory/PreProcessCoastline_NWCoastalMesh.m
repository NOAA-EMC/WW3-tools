
SetPath

MakeCoastalBoundariesOSM.m
lonWest=-128;lonEast=-121;
latSouth=45;latNorth=51;
CoastLineFile = 'GlobalCoastlineOSM.shp'
BuildBoundaryPSLGfunction(CoastLineFile,lonWest,lonEast,latSouth,latNorth)
%[xpi,ypi]=ExtraPointsOfIntrest;
TargetShap=uscl;
DX=.0125;
MakeDistanceToCoastData(DX,TargetShape,[],[])
