
#
# This script partitions the nodes in mshnm into NP parts and writes a script to
# interpolate bathymetry data to the nodes in NP parallel parts
# call:
#
# $ python3 DivideMeshNodes.py 100
#
# to create jobcard: jobcardInterpBathy2Mesh, which is then run: 
#
# $ sbatch jobcardInterpBathy2Mesh
#
# to interpolate the bathymetry in 100 parallel tasks
#
# This script should be run after the bathymetry data sets specified in 
# InterpolateCRM.partscat.py are present
#
"""
Ursa environment

module purge
module use /scratch3/NCEPDEV/climate/Keston.Smith/global-workflowA/sorc/ufs_model.fd/modulefiles
module load ufs_ursa.intel
module load py-scipy/1.14.1
module load py-netcdf4/1.7.1.post2
pip list
"""

import os
import argparse
import time
import numpy as np
import netCDF4 as nc
import math
import sys
from scipy.interpolate import RegularGridInterpolator
#import re
import FiniteElementMeshRoutines as FE
import random

mshnm="HawaiiTest"
mesh="meshes/"+mshnm+".msh"
OutDir=mshnm+".files/"

jobcard="jobcardInterpBathy2Mesh"


def WriteInterpJobscript(fl,N, ComputeNodes):
    with open(fl, 'w') as f:
        #yi[k]=float(values[4])
        f.write("#!/bin/bash \n")
        f.write("#SBATCH --job-name=CRM_interp_masterscript \n")
#        f.write("#SBATCH --ntasks="+str(N)+" \n")
        f.write("#SBATCH --ntasks=1 \n") # ntasks per interpolation
        f.write("#SBATCH --time=08:00:00 \n") 
        f.write("#SBATCH --output=mpi_test_%j.log \n")
        f.write("#SBATCH --error=%j.err \n")
        f.write("#SBATCH --account=marine-cpu \n")
        f.write("#SBATCH --nodes="+str(ComputeNodes)+" \n")
        f.write("#SBATCH --ntasks-per-core=1"+" \n")
        f.write("#SBATCH --array=0-"+str(N-1)+" \n")

        f.write(" \n")

        f.write("module purge \n")
        f.write("module use /scratch3/NCEPDEV/climate/Keston.Smith/global-workflowA/sorc/ufs_model.fd/modulefiles \n")
        f.write("module load ufs_ursa.intel \n")
        f.write("module load py-scipy/1.14.1 \n")
        f.write("module load py-netcdf4/1.7.1.post2 \n")
        f.write("pip list \n")
        f.write("srun python3  InterpolateCRM.partscat.py $SLURM_ARRAY_TASK_ID > InterpJob.$SLURM_ARRAY_TASK_ID.out \n")
        f.write("wait \n")
        f.write("## Run this after the full job array is compleate.\n")
        f.write("## Patch randomly shuffled fields in orgigonal order.\n")
        
        f.write("python3 KnitOutputBackTogether.py "+str(N)+" GMN\n")
        f.write("python3 KnitOutputBackTogether.py "+str(N)+" GMU\n")
        f.write("python3 KnitOutputBackTogether.py "+str(N)+" InvDist\n")
        f.write("python3 KnitOutputBackTogether.py "+str(N)+" NDataPoints\n")
        f.write("python3 KnitOutputBackTogether.py "+str(N)+" ClosestValue\n")
        f.write("python3 KnitOutputBackTogether.py "+str(N)+" LengthScale\n")
        f.write("python3 KnitOutputBackTogether.py "+str(N)+" GMU.std\n")
        f.write("python3 KnitOutputBackTogether.py "+str(N)+" GMN.std\n")
        #f.write("python3 KnitOutputBackTogether.py "+str(N)+" GMM\n")
        #f.write("python3 KnitOutputBackTogether.py "+str(N)+" GM0\n")


# Create the output directory------------------------------------------
try:
    os.mkdir(OutDir)
    print(f"Directory '{OutDir}' created successfully.")
except FileExistsError:
    print(f"Directory '{OutDir}' already exists. Proceeding ...")
except PermissionError:
    print(f"Permission denied: Unable to create '{OutDir}'.")

fl="meshes/"+mshnm+".msh"

xi, yi, ei =FE.loadWW3MeshCoords(fl)

#tmp
lsE=FE.lengthscale(xi, yi, ei)
areaE=FE.ElementArea(xi, yi, ei)
lsN=FE.ComputeNodeLengthScale(lsE, areaE, ei)
np.savetxt("LengthScaleNodes.txt",  lsN, fmt='%.6f', delimiter='\n')

# randomly sort nodes for roughly equal jobs size
# nodes within CRM envelope use far more compute than
# open ocean nodes with fewer points used in interpolation

nn=len(xi)
Nodes = list(range( nn ))
random.shuffle(Nodes) # randomize order of nodes for equitable job load

Nparts=int(sys.argv[1])

nn=xi.shape[0]
NodesPerProc=math.ceil(nn/Nparts)
print(str(NodesPerProc)+ " " + str(nn))
a = range(nn)

NodeList = [] 
for i in range(0, len(a), NodesPerProc):  # Slice list in steps of n
    NodeList.append(a[i:i + NodesPerProc])

print("NodeList")
print(NodeList)
for k in range(Nparts):
    flo=OutDir+'NodeList.'+str(k)+'.txt'
    f=open(flo, 'w')
    f.write("List of nodes for job "+str(k) + "\n")
    List=NodeList[k]
    n=len(List)
    f.write(str(n) + "\n")
    for j in range(n):
        #f.write(str(List[j]) + '\n') # sequential order of nodes, output can be cat > back together
        f.write(str(Nodes[List[j]]) + '\n') # random order of nodes
    f.close

#write jobcard to do interpolation
WriteInterpJobscript(jobcard,Nparts,1)
