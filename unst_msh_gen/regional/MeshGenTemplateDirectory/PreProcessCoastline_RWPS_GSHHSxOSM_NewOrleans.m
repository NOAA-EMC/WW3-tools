% Preprocess inputs to jigsaw to create RWPS style mesh using 
% a coasline based on a blend of OSM and GSHHS with a modification
% near New Orleans to reflect proposed new marine zones for the
% area.  After running this, run RWPSMeshGenScript.GSSHSxOSM.NewOrleans.py 
% to generate jigsaw mesh. Finally run PostProcessGrid_GSHHSxOSM_NewOrleans.m.


SetPath

system('cp /scratch3/NCEPDEV/climate/Keston.Smith/RWPS/Data/JigsawFormatFiles/PSLGboundary1kmP_NewOrleans.GSHHS.msh ./RWPS.PSLG.NewOrleans.GSHHS.msh')
system('cp /scratch3/NCEPDEV/climate/Keston.Smith/RWPS/Data/JigsawFormatFiles/PSLGboundary1kmP_NewOrleans.OSM.msh ./RWPS.PSLG.NewOrleans.OSM.msh')
MergePslgsOSMxGSSHS
%use MakeDistanceToCoastRWPS rather than MakeDistanceToCoast for RWPS to deal with international dateline discontinuity.
MakeDistanceToCoastRWPS



