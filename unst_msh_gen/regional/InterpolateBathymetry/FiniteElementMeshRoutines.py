#
# Functions for reading and writing WW3 meshes, calculating element 
# area, and estimateing mesh length scale.
#
# x and y are length nn node coordinates
# e is a (Number-of-Elements by 3) element array
# bnd is a length nb list of open ocean boundary nodes
#   Note:(nodes are assumed to be indexed from 1,2,...,nn within e
#         and bnd)
#
# loadWW3MeshCoords(fl): read mesh nodes and elements from file
#
# loadWW3Mesh(fl):  read mesh nodes, elements, and open boundary nodes
#                   from WW3 .msh file
#
# WriteWW3Mesh(flo,x,y,z,e,bnd): write WW3 .msh file
#
# lengthscale(x, y, e): compute length scale for each element
#
# ElementArea(x, y, e): compute area for each element
#
# ComputeNodeLengthScale(x,y,e): estimate length scale at each node
#

import numpy as np
import math
#from geopy import distance
import re

# Simple distance approximation in kilometers for calculating element
# area and lengthscale.
deg2km=111.132954
deg2rad=np.pi/180.
def Distance(x0,y0,x1,y1):
    d=np.sqrt( 
        ( deg2km*np.cos(deg2rad*(y0+y1)/2 )*(x1-x0) )**2 + 
        ( deg2km*(y1-y0))**2  ) 
    return d

def loadWW3MeshCoords(fl):
    f=open(fl, 'r')
    header = f.readline() 
    header = f.readline() 
    header = f.readline() 
    header = f.readline() 
    header = f.readline() # number of nodes
    nn = re.findall(r'\d+', header)
    nn=int(nn[0])
    print(nn)
    xi=np.zeros(nn)
    yi=np.zeros(nn)
    k=0
    for i in range(nn):
        A = f.readline()
        values = A.split(" ")
        if len(values)>4:
            xi[k]=values[2]
            yi[k]=values[4]
        else:
            xi[k]=values[1]
            yi[k]=values[2]
        k=k+1
    print("number of nodes read: "+str(k))
    header = f.readline() 
    header = f.readline() 
    header = f.readline() # number of elements
    ne=int(header)
    print("ne=",str(ne))
    ei=np.zeros((ne,3), dtype=int)
    k=0
    for i in range(ne):
        A = f.readline()
        values = A.split(" ")
        #need to factor bnd nodes
        if len(values)>15:
            ei[k,0]=int(values[12])
            ei[k,1]=int(values[14])
            ei[k,2]=int(values[16])
            k=k+1
        elif len(values)>7:
            ei[k,0]=int(values[6])
            ei[k,1]=int(values[7])
            ei[k,2]=int(values[8])
            k=k+1
    print("number of elements read: "+str(k))
    return xi, yi, ei


def loadWW3Mesh(fl):
    f=open(fl, 'r')
    header = f.readline() 
    header = f.readline() 
    header = f.readline() 
    header = f.readline() 
    header = f.readline() # number of nodes
    nn=int(header)
    print("nn = "+str(nn))
    xi=np.zeros(nn)
    yi=np.zeros(nn)
    zi=np.zeros(nn)
    k=0
    for i in range(nn):
        A = f.readline()
        values = A.split(" ")
        if len(values)>4:
            xi[k]=values[2]
            yi[k]=values[4]
            zi[k]=values[6]
        else:
            xi[k]=values[1]
            yi[k]=values[2]
            zi[k]=values[3]
        k=k+1
    print("number of nodes read: "+str(k))
    header = f.readline() 
    header = f.readline() 
    header = f.readline() # number of elements
    ne=int(header)#includes boundary nodes and actual elements
    print("ne="+str(ne)+" -includes boundary nodes")
    nbnd=0
    bnd=[]
    eix=np.zeros((ne,3), dtype=int)
    k=0
    for i in range(ne):
        A = f.readline()
        values = A.split(" ")
        if len(values) == 6:
            if int(values[2])==2:
                bnd.append(int(values[5]))
                nbnd=nbnd+1
        if len(values)>15:
            eix[k,0]=int(values[12])
            eix[k,1]=int(values[14])
            eix[k,2]=int(values[16])
            k=k+1
        elif len(values)>7:
            eix[k,0]=int(values[6])
            eix[k,1]=int(values[7])
            eix[k,2]=int(values[8])
            k=k+1
    ei=eix[range(k),:]
    print("number of open boundary nodes read: "+str(nbnd))
    print("number of elements read: "+str(k))
    return xi, yi, zi, ei, bnd

def WriteWW3Mesh(flo,x,y,z,e,bnd):
    f=open(flo, 'w')
    f.write("$MeshFormat\n")
    f.write("2 0 8\n")
    f.write("$EndMeshFormat\n")
    f.write("$Nodes\n")
    nn=len(x)
    f.write(str(nn)+"\n")
    for k in range(nn):
        f.write(f"{k+1:8d} {x[k]:6f} {y[k]:6f} {z[k]:5f} \n")
    ne=e.shape[0]
    nbnd=len(bnd)
    f.write("$EndNodes\n")
    f.write("$Elements\n")
    f.write(str(ne+nbnd)+"\n")
    for k in range(nbnd):
        f.write(str(k+1)+" 15 2 0 0 " + str(bnd[k])+"\n")
    for k in range(ne):
        f.write(str(nbnd+k+1)+" 2 3 0 "+str(k+1)+" 0 "+str(e[k,0])+" "+str(e[k,1])+" "+str(e[k,2])+"\n")
    f.write("$EndElements\n")
    f.close

# compute length scale for each element e
def lengthscale(x, y, e):
    ne=e.shape[0]
    lengthscaleE=np.zeros(ne)
    for k in range(ne):
        if k % 10000 ==0:
            print(k)
        x1=x[e[k,0]-1]
        y1=y[e[k,0]-1]
        x2=x[e[k,1]-1]
        y2=y[e[k,1]-1]
        x3=x[e[k,2]-1]
        y3=y[e[k,2]-1]
        D3=Distance(x1,y1,x2,y2)
        D1=Distance(x2,y2,x3,y3)
        D2=Distance(x3,y3,x1,y1)
        lengthscaleE[k]=(D1+D2+D3)/3  
        # mean edge length
    return lengthscaleE

# compute area for each element e
def ElementArea(x, y, e):
    ne=e.shape[0]
    AreaE=np.zeros(ne)
    for k in range(ne):
        if k % 10000 ==0:
            print(k)
        x1=x[e[k,0]-1];
        y1=y[e[k,0]-1];
        x2=x[e[k,1]-1];
        y2=y[e[k,1]-1];
        x3=x[e[k,2]-1];
        y3=y[e[k,2]-1];
        D3=Distance(x1,y1,x2,y2)
        D1=Distance(x2,y2,x3,y3)
        D2=Distance(x3,y3,x1,y1)
        S=(D1+D2+D3)/2
        A=np.sqrt(S * (S - D1) * (S - D2) * (S - D3))
        AreaE[k]=A
    return AreaE

# approximate mesh lengthscale at mesh nodes.
# lengthscale is the weighted average (by area) of surrounding elements
def ComputeNodeLengthScale(x,y,e):
    nn=np.max(e[:])
    print("ComputeNodeLengthScale nn="+str(nn))
    areaE=ElementArea(x, y, e)
    LengthScaleE=lengthscale(x, y, e)
    
    ne=e.shape[0]
    areaT=np.zeros(nn)
    lsN=np.zeros(nn)
    for k in range(ne):
        for j in range(3):
            i=e[k,j]-1
            lsN[i]=lsN[i]+LengthScaleE[k]*areaE[k]
            areaT[i]=areaT[i]+areaE[k]
    for k in range(nn):
        lsN[k]=lsN[k]/areaT[k]
    return lsN

