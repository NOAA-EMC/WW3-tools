
%Preprocess jigsaw inputs for building a RWPS type mesh using the GSHHS coastline

SetPath
MakeCoastalBoundariesGSHHS 
lonWest=129.91;lonEast=10.71;latSouth=-30.42;latNorth=79.99;
CoastLineFile = 'GlobalCoastlineGSHHS.shp'
BuildBoundaryPSLGfunction(CoastLineFile,lonWest,lonEast,latSouth,latNorth)

% use MakeDistanceToCoastRWPS, rather than MakeDistanceToCoastData, for
% RWPS to deal with  international dateline discontinuity.
MakeDistanceToCoastRWPS 



