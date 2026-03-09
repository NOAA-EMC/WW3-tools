% This file is used to set global variable values used in the regional
% mesh generation routines.  It also adds matlab libraries used by the
% routines to the matlab path.  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Define global variables used in regional mesh generation %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global GlobalTopoFile TargetCoastlineFile GlobalCoastlineFile 
global GlobalCoastlineFileGSHHS GlobalCoastlineFileOSM PSLGfile 
global LakeDir NewOrleansCoastDir
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% STATIC FILE PATH VARIABLES %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These files are downloaded by script:
% unst_msh_gen/regional/MeshGenTemplateDirectory/DownloadShapeFiles.sh
% and should not need to be changed to run demonstration cases.

% Global coverage "low" resolution bathymetry file
GlobalTopoFile='../Data/Bathymetry/RTopo_2_0_4_GEBCO_v2024_60sec_pixel.nc'

% US coastline file or shapefile defining coastal points where high
% resolution in the mesh is desired 
TargetCoastlineFile='../Data/us_coastline/tl_2023_us_coastline.shp'

%Global coverage coastline file from Global Self-consistent, 
%Hierarchical, High-resolution Geography Database (GSHHG) 
GlobalCoastlineFileGSHHS='../Data/GSHHS_shp/f/GSHHS_f_L1.shp'

%Global coverage coastline file from Open Street Map (OSM)
GlobalCoastlineFileOSM='../Data/openstreetmap_land/land-polygons-complete-4326/land_polygons.shp'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%!!!!!!!!! EDIT BELOW HERE TO SWITCH BETWEEN DEMO CASES !!!!!!!!!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% DEFINE CASE SPECIFIC FILE PATHS HERE %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Filename "PSLGfile" points to a file that is created by calling 
% PreProcessCoastline_X.m and defines the boundary geometry of the 
% mesh as a Planar straight Line Graph (PSLG). This file is output 
% from function BuildBoundaryPSLGfunction and must match file 
% specification in python script MeshGenScript.X.py. The file is 
% also an input for function MakeDistanceToCoast.
%
% Filename "GlobalCoastlineFile" specifies which global coastline file
% is used to designate land-ocean boundaries in the creation of the 
% mesh boundary. Shapefiles for the local region could also be used.

% For NWcoastal demonstration mesh uncomment the following lines:
GlobalCoastlineFile=GlobalCoastlineFileOSM
PSLGfile='NWcoastal.PSLG.msh'

% For RWPS type mesh with GSHHS coastline set:
%PSLGfile='GlobalCoastlineGSHHS.PSLG.msh'
%GlobalCoastlineFile=GlobalCoastlineFileGSHHS

% For RWPS mesh with modified New Orleans coastline set:
%PSLGfile='RWPS.GSHHSxOSM.NewOrleans.PSLG.msh'
% this case does not require specification of GlobalCoastlineFile
% as it builds the boundary from files located in directory 
% "NewOrleansCoastDir"(see below)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%!!!!!!!!! EDIT ABOVE HERE TO SWITCH BETWEEN DEMO CASES !!!!!!!!!!!!!!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% RELEVENT TO RWPS MESH ONLY %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These variables are only used for the RWPS domain. The "LakeDir" 
% contains existing meshes representing inland lakes that are to be 
% included in the RWPS mesh. The "NewOrleansCoastDir" is the path to a 
% directory containing .msh files with a modified coastline around New
% Orleans defined by proposed New Marine Zones for the region.

% On Ursa:
LakeDir='/scratch3/NCEPDEV/climate/Keston.Smith/RWPS/RWPSLakes/' %Ursa
NewOrleansCoastDir='/scratch3/NCEPDEV/climate/Keston.Smith/RWPS/Data/JigsawFormatFiles/' %Ursa

% On Orion:
%LakeDir='/work2/noaa/marine/keston/DATA/NWPSLakes/' % Orion
%NewOrleansCoastDir='/work2/noaa/marine/keston/DATA/RWPSmshFiles/' % Orion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set matlab path to include matlab libraries used here. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath ../matlab
addpath ../matlab/graphics
addpath ../matlab/jigsaw-matlab
addpath ../matlab/inpoly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
