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

3- clone the repo
        $git clone https://github.com/NOAA-EMC/WW3-tools
	$cd WW3-tools/unst_msh_gen

3- get the DEM and make sure it is in the WW3-tools/unst_msh_gen directory
	$wget https://github.com/dengwirda/dem/releases/download/v0.1.1/RTopo_2_0_4_GEBCO_v2023_60sec_pixel.zip
4- unzip the DEMi 
	unzip *.zip
 
## Usage
5- run the script inside of WW3-tools/unst_msh_gen:
	$python3 ocn_ww3.py --black_sea [option]
		
		option=  3: default which will have the Black Sea and the connections.
			 2: will have the Balck sea as a seperate basin.
			 1: will exclude the Black sea

NOTE: the output will be gmsh format which will be used by WW3.

NOTE: for different resolution (uniform) in km, you should change the following:
	opts.hfun_hmax
	hmax
	hshr
	hmin

NOTE: The output mesh will have -180:180 longitude, you can convert this by unisg ShiftMesh.py script, to 0:360 longitude.
	input_file_path: your jigsaw mesh in gmsh format with -18:180 long
	output_file_path: shifted mesh in gmsh format with 0:360 long


## Contributing
This is ongoing effort with the great help of Darren Engwirda, JIGSAW developer.
