#!/bin/bash 
#############################################################################################################

# Script to download global and regional bathymetry files:
# You will need ~ 20 GB free space

echo "Downloading bathymetry data"

#############################################################################################################
######################## Download main Coastal Relief Moesl (CRM) files #####################################
#############################################################################################################
# see: https://www.ncei.noaa.gov/products/coastal-relief-model
#############################################################################################################
#Northeast Atlantic
wget --output-document crm_vol1_2023.nc https://www.ngdc.noaa.gov/thredds/fileServer/crm/cudem/crm_vol1_2023.nc
#Southeast Atlantic
wget --output-document crm_vol2_2023.nc https://www.ngdc.noaa.gov/thredds/fileServer/crm/cudem/crm_vol2_2023.nc
# Florida and East Gulf of America
wget --output-document crm_vol3_2023.nc https://www.ngdc.noaa.gov/thredds/fileServer/crm/cudem/crm_vol3_2023.nc
# Central GOA
wget --output-document crm_vol4_2023.nc https://www.ngdc.noaa.gov/thredds/fileServer/crm/cudem/crm_vol4_2023.nc
# Western GOA
wget --output-document crm_vol5_2023.nc https://www.ngdc.noaa.gov/thredds/fileServer/crm/cudem/crm_vol5_2023.nc
#Central Pacific
wget --output-document crm_vol7_2023.nc https://www.ngdc.noaa.gov/thredds/fileServer/crm/cudem/crm_vol7_2025.nc
#Northwest Pacific
wget --output-document crm_vol8_2023.nc https://www.ngdc.noaa.gov/thredds/fileServer/crm/cudem/crm_vol8_2025.nc
#Puerto Rico
wget --output-document crm_vol9_2023.nc https://www.ngdc.noaa.gov/thredds/fileServer/crm/cudem/crm_vol9_2023.nc
#Hawaii
wget --output-document crm_vol10_2023.nc https://www.ngdc.noaa.gov/thredds/fileServer/crm/cudem/crm_vol10_2023.nc
#Other CRM files
# Southern California
wget --output-document crm_vol6.nc https://www.ngdc.noaa.gov/thredds/fileServer/crm/crm_vol6.nc
# Southern Alaska
wget --output-document crm_southak.nc https://www.ngdc.noaa.gov/thredds/fileServer/crm/crm_southak.nc
#############################################################################################################

echo "Done Coastal Relief Model (CRM) Download"

#############################################################################################################
######################## Download regional files in Pacific #################################################
#############################################################################################################
# see: https://pae-paha.pacioos.hawaii.edu/thredds/bathymetry.html
#############################################################################################################
# American Samoa
wget --output-document ngdc_bathy_90m_amsamoa.nc https://pae-paha.pacioos.hawaii.edu/thredds/ncss/ngdc_bathy_90m_amsamoa?var=elev&disableLLSubset=on&disableProjSubset=on&horizStride=1&addLatLon=true
# Mariana islands
wget --output-document ngdc_bathy_180m_mariana.nc https://pae-paha.pacioos.hawaii.edu/thredds/ncss/ngdc_bathy_180m_mariana?var=elev&disableLLSubset=on&disableProjSubset=on&horizStride=1&addLatLon=true
#North Western Hawiian islands
wget --output-document hurl_bathy_60m_nwhi.nc https://pae-paha.pacioos.hawaii.edu/thredds/ncss/hurl_bathy_60m_nwhi?var=elev&disableLLSubset=on&disableProjSubset=on&horizStride=1&addLatLon=true
#############################################################################################################

echo "Done regional files in Pacific Download"

#############################################################################################################
######################## Download small scale bathymetry for US Minor outlying islands in Pacific ###########
#############################################################################################################
# see: https://pae-paha.pacioos.hawaii.edu/thredds/bathymetry.html
#############################################################################################################
wget --output-document ngdc_bathy_10m_wake.nc https://pae-paha.pacioos.hawaii.edu/thredds/ncss/ngdc_bathy_10m_wake?var=elev&disableLLSubset=on&disableProjSubset=on&horizStride=1&addLatLon=true

wget --output-document pibhmc_bathy_20m_jarvis.nc https://pae-paha.pacioos.hawaii.edu/thredds/ncss/pibhmc_bathy_20m_jarvis?var=elev&disableLLSubset=on&disableProjSubset=on&horizStride=1&addLatLon=true

wget --output-document pibhmc_bathy_20m_johnston.nc https://pae-paha.pacioos.hawaii.edu/thredds/ncss/pibhmc_bathy_20m_johnston?var=elev&disableLLSubset=on&disableProjSubset=on&horizStride=1&addLatLon=true

wget --output-document pibhmc_bathy_20m_kingman.nc https://pae-paha.pacioos.hawaii.edu/thredds/ncss/pibhmc_bathy_20m_kingman?var=elev&disableLLSubset=on&disableProjSubset=on&horizStride=1&addLatLon=true

wget --output-document pibhmc_bathy_40m_baker.nc https://pae-paha.pacioos.hawaii.edu/thredds/ncss/pibhmc_bathy_40m_baker?var=elev&disableLLSubset=on&disableProjSubset=on&horizStride=1&addLatLon=true

wget --output-document pibhmc_bathy_40m_howland.nc https://pae-paha.pacioos.hawaii.edu/thredds/ncss/pibhmc_bathy_40m_howland?var=elev&disableLLSubset=on&disableProjSubset=on&horizStride=1&addLatLon=true

wget --output-document pibhmc_bathy_40m_palmyra.nc https://pae-paha.pacioos.hawaii.edu/thredds/ncss/pibhmc_bathy_40m_palmyra?var=elev&disableLLSubset=on&disableProjSubset=on&horizStride=1&addLatLon=true

wget --output-document pibhmc_bathy_40m_rose.nc https://pae-paha.pacioos.hawaii.edu/thredds/ncss/pibhmc_bathy_40m_rose?var=elev&disableLLSubset=on&disableProjSubset=on&horizStride=1&addLatLon=true

wget --output-document pibhmc_bathy_40m_vailuluu.nc https://pae-paha.pacioos.hawaii.edu/thredds/ncss/pibhmc_bathy_40m_vailuluu?var=elev&disableLLSubset=on&disableProjSubset=on&horizStride=1&addLatLon=true

wget --output-document pibhmc_bathy_40m_swains.nc https://pae-paha.pacioos.hawaii.edu/thredds/ncss/pibhmc_bathy_40m_swains?var=elev&disableLLSubset=on&disableProjSubset=on&horizStride=1&addLatLon=true

wget --output-document pibhmc_bathy_5m_palmyra.nc https://pae-paha.pacioos.hawaii.edu/thredds/ncss/pibhmc_bathy_5m_palmyra?var=elev&disableLLSubset=on&disableProjSubset=on&horizStride=1&addLatLon=true

wget --output-document sopac_bathy_50m_majuro_reef.nc https://pae-paha.pacioos.hawaii.edu/thredds/ncss/sopac_bathy_50m_majuro_reef?var=elev&disableLLSubset=on&disableProjSubset=on&horizStride=1&addLatLon=true
#############################################################################################################

echo "Done small scale bathymetry for US Minor outlying islands in Pacific download"

#############################################################################################################
######################## Global SDB [~0m to ~30m] depth from copernicus #####################################
#############################################################################################################
# see: https://data.marine.copernicus.eu/product/BATHYMETRY_GLO_PHY_COASTAL_L4_MY_016_001/files?subdataset=cmems_obs-sdb_glo_phy_comp_my_100m-l4-s2_static_202511
#############################################################################################################
wget --output-document cmems_obs-sdb_glo_phy-comp_my-oa-100m-l4-s2_static.nc https://s3.waw3-1.cloudferro.com/mdl-native-17/native/BATHYMETRY_GLO_PHY_COASTAL_L4_MY_016_001/cmems_obs-sdb_glo_phy_wk_my_100m-l4-s2_static_202511/cmems_obs-sdb_glo_phy-wk_my-oa-100m-l4-s2_static.nc
#############################################################################################################

echo "Done Global SDB download"

#############################################################################################################
######################## 60 second Global Bathymetry#########################################################
wget https://github.com/dengwirda/dem/releases/download/v0.1.1/RTopo_2_0_4_GEBCO_v2024_60sec_pixel.zip
#############################################################################################################

echo "Done Global Bathymetry download"
