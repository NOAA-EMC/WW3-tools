
#
# This script takes the files generated in parallel from  InterpolateCRM.partscat.py
# It should be called at the end of 
# $ sbatch jobcardInterpBathy2Mesh
# when all independent interpolation problems are finished
# 

import numpy as np
import FiniteElementMeshRoutines as FE
import sys

mshnm="HawaiiTest"
mesh="meshes/"+mshnm+".msh"
OutDir=mshnm+".files/"

Nparts=int(sys.argv[1])
Prefix=sys.argv[2]
xi, yi, ei = FE.loadWW3MeshCoords(mesh)
nn=len(xi)
NodeList=[]
si=np.zeros(nn)
for n in range(Nparts):
    fli=OutDir+'NodeList.'+str(n)+'.txt'
    fpi=open(fli, 'r')
    fls=OutDir+Prefix+"."+str(n)+".txt"
    fps=open(fls, 'r')
    print("reading nodelist from: "+fli+", reading values from: "+fls)
    header = fpi.readline()
    nn0=int(fpi.readline()) 
    for k in range(nn0):
        j = int(fpi.readline())
        value=float(fps.readline())
        si[j] = value
fps.close
fpi.close

#write full field file with interpolated values
flout=mshnm+"."+Prefix+".txt"
np.savetxt(flout , si, fmt='%.6f',delimiter='\n')

