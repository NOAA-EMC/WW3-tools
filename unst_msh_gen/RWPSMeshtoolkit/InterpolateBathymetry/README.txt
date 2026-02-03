These are subroutines to interpolate bathymetry data to a unstructured
meshes with triangular elements.  The interpolation uses Gauss Markov 
smoothing with variable lengthscale proportional to the mesh lengthscale.
The interpolation is computed in parallel across a user determined 
number of mesh divisions.  The program uses an arbitrary number of 
bathymetry data files.  It was designed to use a low resolution file 
with global coverage as well as local high resolution gridded bathymetry,
and global scattered high resolution bathymetry data, such as global
Satellite Derived Bathymetry (SDB). Bathymetry data sets need to be 
downloaded and preprocessed to use a common format. 

To interpolate bathymetry to a mesh:

    (0) Download bathymetry data to be used (notes are in  
        MakeCommonNetCDFFileFormat.py)
       
    (1) python3 MakeCommonNetCDFFileFormat.py

    (2) python3 SmoothSubsampleNetCDFFiles.py

    (3) python3 SmoothSubsampleCMEMS.py

    (4) python3 DivideMeshNodes.py NP
            NP=number of interpolation jobs to run in parallel

    (5) sbatch jobcardInterpBathy2Mesh
            Runs interpolation job in parallel in slurm environmant
    
    (6) python3 AddBathyToMesh.py HawaiiTest.GMU.txt
            creates HawaiiTest.GMU.txt.WW3.msh, a WW3 .msh mesh format

to change the mesh interpolated to :

$sed -i 's/HawaiiTest/MyMesh/g' *.py

and copy MyMesh.msh file to the 
WW3-tools/unst_msh_gen/RWPSMeshtoolkit/InterpolateBathymetry/meshes/ 
directory


AddBathyToMesh.py
FiniteElementMeshRoutines.py
GaussMarkov.py
InterpolateCRM.partscat.py
KnitOutputBackTogether.py
MakeCommonNetCDFFileFormat.py
SmoothSubsampleCMEMS.py
SmoothSubsampleNetCDFFiles.py
