#!/bin/bash 

# Script to downoad files used in mesh generation for RWPS 
echo "Downloading coastline shapefiles for boundary definition"

#############################################################################################################
######################## Download Coastline Shape Files ####################################################
#############################################################################################################

#US coastline data
echo "Downloading US coastline"
wget --output-document tl_2023_us_coastline.zip https://www2.census.gov/geo/tiger/TIGER2023/COASTLINE/tl_2023_us_coastline.zip

echo "Downloading Global Self-consistent, Hierarchical, High-resolution Shorelines (GSHHS)"
#GSHHG global coastline data
wget --output-document gshhg-shp-2.3.7.zip http://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.7.zip

echo "Downloading OpenStreetMap shoreline (OSM)"
#OSM global coastline data
wget --output-document land-polygons-complete-4326.zip https://osmdata.openstreetmap.de/download/land-polygons-complete-4326.zip

#############################################################################################################
######################## 60 second Global Bathymetry#########################################################

echo "Downloading global bathymetry"
wget https://github.com/dengwirda/dem/releases/download/v0.1.1/RTopo_2_0_4_GEBCO_v2024_60sec_pixel.zip

#############################################################################################################
