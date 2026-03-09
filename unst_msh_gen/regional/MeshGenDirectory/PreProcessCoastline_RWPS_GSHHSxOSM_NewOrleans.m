% Preprocess inputs to jigsaw to create RWPS style mesh using 
% a coasline based on a blend of OSM and GSHHS with a modification
% near New Orleans to reflect proposed new marine zones for the
% area.  After running this, run RWPSMeshGenScript.GSHHSxOSM.NewOrleans.py 
% to generate jigsaw mesh. Finally run PostProcessGrid_GSHHSxOSM_NewOrleans.m.

SetPath

system(['cp ',NewOrleansCoastDir,'PSLGboundary1kmP_NewOrleans.GSHHS.msh ./RWPS.PSLG.NewOrleans.GSHHS.msh']);
system(['cp ',NewOrleansCoastDir,'PSLGboundary1kmP_NewOrleans.OSM.msh ./RWPS.PSLG.NewOrleans.OSM.msh']);

% Substitute GSHHS with OSM coastlines where problems are location 
% errors have been found in GSHHS, for example American Samoa.
MergePslgsOSMxGSHHS 

% use MakeDistanceToCoastRWPS, rather than MakeDistanceToCoastData, for
% RWPS to deal with  international dateline discontinuity.
MakeDistanceToCoastRWPS



