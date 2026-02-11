
addpath ../matlab
addpath ../matlab/jigsaw-matlab
addpath ../matlab/inpoly

%addpath /scratch3/NCEPDEV/climate/Keston.Smith/MeshGenMatlabLibs/jigsaw-matlab
%addpath /scratch3/NCEPDEV/climate/Keston.Smith/MeshGenMatlabLibs/InsidePoly


% set global variables defining paths to different data sets
global GlobalTopoFile uscl gcfl gcflGSHHS gcflOSM PSLGfile

% Global coverage "low" resolution bathymetry file
GlobalTopoFile='../Data/Bathymetry/RTopo_2_0_4_GEBCO_v2024_60sec_pixel.nc'

% US coastline file or shapefile defining coastal points for high resolution 
uscl='../Data/us_coastline/tl_2023_us_coastline.shp'

%Global coverage coastline file from ()
gcflGSHHS='../Data/GSHHS_shp/f/GSHHS_f_L1.shp'

%Global coverage coastline file from OpenStreetMap
gcflOSM='../Data/openstreetmap_land/land-polygons-complete-4326/land_polygons.shp'

%which Global coverage coastline file  is to be used
gcfl=gcflGSHHS

%Output file from BuildBoundaryPSLGwGSHHS and input to MakeDfunGSHHS
PSLGfile='GlobalCoastlineGSHHS.PSLG.msh'

