
SetPath

MakeCoastalBoundariesOSM
lonWest=-130;lonEast=-121;
latSouth=45;latNorth=51;
CoastLineFile = 'GlobalCoastlineOSM.shp'
BuildBoundaryPSLGfunction(CoastLineFile,lonWest,lonEast,latSouth,latNorth,PSLGfile)
%[xpi,ypi]=ExtraPointsOfIntrest;
TargetShape=uscl;% US coastline as target for distance
DX=.0125;
MakeDistanceToCoastData(DX,uscl,[],[])
