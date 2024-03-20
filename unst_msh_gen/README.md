# Project Name
Unstructured mesh generation for WW3 using JIGSSAW.
## Description
Mesh generation script capable of creating unstructured meshes for WW3 global modeling. The tool leverages JIGSAWPY (https://github.com/dengwirda/jigsaw-python) for efficient triangulation.
Main changes include:
Implementation of the ocn_ww3.py to create uniform unstructured mesh for global WW3 model.
This tool is under active development, with future work focused on variable unstructured mesh generation.

## Installation

1- Install jigsawpy (https://github.com/dengwirda/jigsaw-python)
2- you need following packages:
	numpy
        scipy
	packaging
	netcdf4
3- get the DEM:
	$wget https://github.com/dengwirda/dem/releases/download/v0.1.1/RTopo_2_0_4_GEBCO_v2023_60sec_pixel.zip
4- unzip the DEM 
	unzip *.zip
 


## Usage
1- clone the repo
	$mkdir meshgen
	$cd meshgen
	$git clone https://github.com/AliS-Noaa/WW3-tools/tree/add-unst-mesh-gen
2- get the DEM:
        $wget https://github.com/dengwirda/dem/releases/download/v0.1.1/RTopo_2_0_4_GEBCO_v2023_60sec_pixel.zip
        $unzip *.zip
3- run the script:
	$python3 ocn_ww3.p --black_sea [option]
		
		option=  3: default which will have the Black Sea and the connections.
			 2: will have the Balck sea as a seperate basin.
			 1: will exclude the Black sea

NOTE: the output will be gmsh format which will be used by WW3.

## Contributing
This is ongoing effort with the great help of Darren Engwirda, JIGSAW developer.
