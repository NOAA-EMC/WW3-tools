

SetPath

%MakeCoastalBoundariesOSM
MakeCoastalBoundariesGSHHSX
%BuildBoundaryPSLGwGSHHS
%MakeDfunGSHHS
SetPath
lonWest=129.91;lonEast=10.71;latSouth=-30.42;latNorth=79.99;
CoastLineFile = 'GlobalCoastlineGSHHS.shp'
BuildBoundaryPSLGfunction(CoastLineFile,lonWest,lonEast,latSouth,latNorth)

%TargetShap=uscl; %Shapefile with points of interest
%DX=.1; %lat/lon grid resolution for jigsaw inputs
%[xpi,ypi]=ExtraPointsOfIntrest; %Define points of interest for distance function not in TargetShape file
%MakeDistanceToCoastData(DX,TargetShape,xpi,ypi) %use for local meshes that don't span the international dateline 

%use this for RWPS to deal with  international dateline discontinuity.
MakeDistanceToCoastRWPS 



