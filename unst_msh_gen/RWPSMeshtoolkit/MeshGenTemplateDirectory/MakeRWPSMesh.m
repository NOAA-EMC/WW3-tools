
SetPath
lonWest=129.91;lonEast=10.71;latSouth=-30.42;latNorth=79.99;
CoastLineFile = 'GlobalCoastlineGSHHS.shp'
BuildBoundaryPSLGfunctionLC(CoastLineFile,lonWest,lonEast,latSouth,latNorth)
[xpi,ypi]=ExtraPointsOfIntrest;
TargetShap=uscl;
DX=.1;
MakeDistanceToCoastData(DX,TargetShape,xpi,ypi)
