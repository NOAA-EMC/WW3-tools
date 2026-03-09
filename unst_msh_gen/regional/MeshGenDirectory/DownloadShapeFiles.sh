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
#wget --output-document gshhg-shp-2.3.7.zip http://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.7.zip
wget --output-document gshhg-shp-2.3.7.zip https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/gshhg-shp-2.3.7.zip

echo "Downloading OpenStreetMap shoreline (OSM)"
#OSM global coastline data
wget --output-document land-polygons-complete-4326.zip https://osmdata.openstreetmap.de/download/land-polygons-complete-4326.zip

#############################################################################################################
######################## 60 second Global Bathymetry#########################################################

echo "Downloading global bathymetry"
wget https://github.com/dengwirda/dem/releases/download/v0.1.1/RTopo_2_0_4_GEBCO_v2024_60sec_pixel.zip

#############################################################################################################

mkdir ../Data
mkdir ../Data/us_coastline/
unzip tl_2023_us_coastline.zip
mv tl_2023_us_coastline.*  ../Data/us_coastline/
echo "US coastline data moved to ../Data/us_coastline/"

unzip gshhg-shp-2.3.7.zip
mv GSHHS_shp/ ../Data/
mv WDBII_shp ../Data/
mv gshhg-shp-2.3.7.zip ../Data
echo "GSHHS global coastline data moved to ../Data/GSHHS_shp/"

unzip land-polygons-complete-4326.zip
mkdir ../Data/openstreetmap_land
mv land-polygons-complete-* ../Data/openstreetmap_land/
echo "OSM global coastline data moved to ../Data/openstreetmap_land/"

unzip RTopo_2_0_4_GEBCO_v2024_60sec_pixel.zip
mkdir ../Data/Bathymetry
mv RTopo_2_0_4_GEBCO_v* ../Data/Bathymetry/
echo "Global bathymetry data moved to ../Data/Bathymetry/"
