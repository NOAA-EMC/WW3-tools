######################################################################
# This is the main interpolation engine.  This script reads in a list 
# of nodes created by DivideMeshNodes.py and interpolates bathymetry 
# from a list of bathymetry data files to the nodes in the assigned
# list.  The intent is to run this in NPart parallel jobs for a 
# partition of mesh nodes into NPart sets.
#
# The user is responsible for specifing the bathymetry data sets
# in list "flnms" as well as the type of data each of these files 
# in list "filetype". filetype[k] = 0 gridded local data, 1 gridded 
# global data and 2 for scattered data. 
#
# The user should define there observational error model and bacground
# process statistical model. 
#
import os
import argparse
import time
import numpy as np
import netCDF4 as nc

#import jigsawpy
import sys
from scipy.interpolate import RegularGridInterpolator
import math

#from geopy import distance
import itertools
#from geopy.distance import haversine
#from pykrige.ok import OrdinaryKriging

import FiniteElementMeshRoutines as FE
import GaussMarkov as GM

# Name of mesh to interpolate bathymetry to
mshnm="HawaiiTest"

mesh="meshes/"+mshnm+".msh"
OutDir=mshnm+".files/"

# TextOutput = True # prints lots of information durring interpolation (may cause slowdown)
TextOutput = False

#  directory containg all of the bahymetry data files specified in flnms
BathyDir="/scratch3/NCEPDEV/climate/Keston.Smith/CoastalReliefModel/"

########################################################
# processed netcdf bathymetry files with variables:
# lon, lat and z. Specify corresponding file type, gridded or scattered,
# in list "filetype" below.
########################################################
flnms=[
    "cmems_obs-sdb_glo_phy-comp_my-oa-100m-l4-s2_static.PointValues500m.nc",
    "crm_vol1_2023.nc.S250m.VB.nc",
    "crm_vol2_2023.nc.S250m.VB.nc",
    "crm_vol3_2023.nc.S250m.VB.nc",
    "crm_vol4_2023.nc.S250m.VB.nc",
    "crm_vol5_2023.nc.S250m.VB.nc",
    "crm_vol7_2024.nc.S250m.VB.nc",
    "crm_vol9_2023.nc.S250m.VB.nc",
    "crm_vol10_2023.nc.S250m.VB.nc",
    "crm_vol6_2023.nc.S250m.VB.nc",
    "crm_vol8_2023.nc.S250m.VB.nc",
    "crm_southak.CRMformat.nc",
    "hurl_bathy_60m_nwhi.CRMformat.nc.S250m.VB.nc",
    "PIX/ngdc_bathy_10m_wake.CRM.nc",
    "PIX/pibhmc_bathy_20m_jarvis.CRM.nc",
    "PIX/pibhmc_bathy_20m_johnston.CRM.nc",
    "PIX/pibhmc_bathy_20m_kingman.CRM.nc",
    "PIX/pibhmc_bathy_40m_baker.CRM.nc",
    "PIX/pibhmc_bathy_40m_howland.CRM.nc",
    "PIX/pibhmc_bathy_40m_palmyra.CRM.nc",
    "PIX/pibhmc_bathy_40m_rose.CRM.nc",
    "PIX/pibhmc_bathy_40m_swains.CRM.nc",
    "PIX/pibhmc_bathy_40m_vailuluu.CRM.nc",
    "PIX/pibhmc_bathy_5m_palmyra.CRM.nc",
    "PIX/sopac_bathy_50m_majuro_reef.CRM.nc",
    "PIX/ngdc_bathy_180m_mariana.CRM.nc",
    "ngdc_bathy_90m_amsamoa.crm.nc",
    "RTopo_2_0_4_GEBCO_v2023_60sec_pixel.CRMformat.nc"
    ]
########################################################

########################################################
# Set filetype corresponding to files in flnms 
# 0 : for gridded local data and 
# 1 : for gridded global data and 
# 2 : for scattered global data
########################################################
nf=len(flnms)
filetype=np.zeros(nf, dtype=int)
filetype[0]=2
filetype[nf-1]=1
########################################################

Zmin=-11000. #deepest legit ocean depth value in case of mask fail 
# nescesarry for data sets "crm_vol6_2023.nc.S250m.nc" and "crm_vol8_2023.nc.S250m.nc"
Zmax=0. # maximum value to include in interpolation - zero out land values or ignore land 
Dmin=5./1000.# 5 meter minimum distance between observation points in km, prevent singularity
lambdaLL=.025 # set deg lat, lon search width for overlapping regions
#use Approximately NxTarget**2 points per data set, 
#each linear system is approx 2*(NxTarget**2)
NxTarget=20 # N=20 ~ 10 grid points in each cardinal direction used in interpolating 
NpointsMax=1000
xlist=[]
ylist=[]
n=0
nxl=np.zeros(len(flnms))
nyl=np.zeros(len(flnms))
xmax=np.zeros(len(flnms))
xmin=np.zeros(len(flnms))
ymax=np.zeros(len(flnms))
ymin=np.zeros(len(flnms))
for fl in flnms:
#    print(f"'{fl}' has a length of {len(fl)}")
#    if n>0:
    data = nc.Dataset(BathyDir+fl,"r")
    x=np.array(data["lon"][:])
    x=x%360
    y=np.array(data["lat"][:])
    xlist.append(list(x))
    ylist.append(list(y))
    nxl[n]=len(x)
    nyl[n]=len(y)
        #setup quick lookup table for file exlusion
    xmax[n]=np.max(x)+lambdaLL*2
    xmin[n]=np.min(x)-lambdaLL*2
    ymax[n]=np.max(y)+lambdaLL*2
    ymin[n]=np.min(y)-lambdaLL*2
    n=n+1

def ZeroPadIntStr(N,K):
    ZPNs=str(N).zfill(K)
    return ZPNs


xi, yi, ei = FE.loadWW3MeshCoords(mesh)

lsN=FE.ComputeNodeLengthScale(xi,yi,ei)

print(np.min(lsN))
print(np.mean(lsN))
print(np.max(lsN))

N=int(sys.argv[1])

fln=OutDir+'NodeList.'+str(N)+'.txt'
floutID=OutDir+'InvDist.'+str(N)+'.txt'
floutGMM=OutDir+'GMM.'+str(N)+'.txt'
floutGMN=OutDir+'GMN.'+str(N)+'.txt'
floutGM0=OutDir+'GM0.'+str(N)+'.txt'
floutGMU=OutDir+'GMU.'+str(N)+'.txt'
floutClosest=OutDir+'ClosestValue.'+str(N)+'.txt'
floutNpts=OutDir+'NDataPoints.'+str(N)+'.txt'
floutLLS=OutDir+'LengthScale.'+str(N)+'.txt'
floutGMUerr= OutDir+'GMU.std.'+str(N)+'.txt'
floutGMNerr= OutDir+'GMN.std.'+str(N)+'.txt'
#read in nodes to interpolate to
f=open(fln,"r")
header = f.readline()
nn =int( f.readline())

#make local node set to interpolate to
xil=np.zeros(nn)
yil=np.zeros(nn)
LSl=np.zeros(nn)
for k in range(nn):
    j = int(f.readline())
    xil[k]=xi[j]
    yil[k]=yi[j]
    LSl[k]=lsN[j]

nn=len(xil)

# set search width for each data set based on expected
# number of points to interpolate
SearchWidth=np.zeros(nf)
for j in range(nf):
    x=np.array(xlist[j][:])
    dx=np.abs(x[1]-x[0])
    SearchWidth[j]=dx*float(NxTarget)/2.
    print("search width(deg Lat lon) for file: "+flnms[j]+" = ",str(SearchWidth[j]))

SearchWidth[0]=.05
SearchWidth[nf-1]=.75*SearchWidth[nf-1] # smaller number for global bathy set
    
# Initialize estimate fields and posterioir estimates
NumPoints=np.zeros(nn)
ziID=np.zeros(nn)
ziGMU=np.zeros(nn)
stdiGMU=np.zeros(nn)
ziGMN=np.zeros(nn)
stdiGMN=np.zeros(nn)
#ziGMM=np.zeros(nn)
#ziGM0=np.zeros(nn)
ziClosest=np.zeros(nn)
LocalLengthScale=np.zeros(nn)

######################################################################
# For each scattered data set, read in the scattered data here to
# one dimensional arrays x0, y0, and z0.  This should be looped to
# accumulate all of the scatted datasets together. The only scattered
# data set used here is
# "cmems_obs-sdb_glo_phy-comp_my-oa-100m-l4-s2_static.PointValues500m.nc",
# so this is simply loaded here.
######################################################################
fl=flnms[0]                
data = nc.Dataset(BathyDir+fl,"r")
x0=data["lon"][:]
x0=x0%360
y0=data["lat"][:]
z0=data["z"][:]
######################################################################

t0 = time.time()
for n in range(nn):
    xp=xil[n]%360
    yp=yil[n]
# Set local length scale to use in Gauss Markov smoothing at node n
#    LSp=2.*LSl[n] #probably better choice for spherical covarianve function
    LSp=1.*LSl[n] #choice for exponential covarianve function
    LocalLengthScale[n]=LSp 
    xs=[]
    ys=[]
    zs=[]
    si=0
    for j in range(nf):
        if filetype[j] <2 : #gridded data
            if all([ xp < xmax[j],xp > xmin[j],yp < ymax[j],yp > ymin[j]]):
                x=np.array(xlist[j][:])
                y=np.array(ylist[j][:])
                fl=flnms[j]
                if TextOutput:
                    print(fl)
                jx=np.array(np.where( np.abs(xp-x) < SearchWidth[j] ))
                jx=list(itertools.chain.from_iterable(jx))
                jy=np.array(np.where( np.abs(yp-y) < SearchWidth[j] ))
                jy=list(itertools.chain.from_iterable(jy))
                data = nc.Dataset(BathyDir+fl,"r")
                for kx in jx:
                    for ky in jy:
                            if kx and ky: # not empty 
                                z=data["z"][ky,kx]
                                if not z.mask:
                                    IncludePoint=True
                                    if len(xs)>0: #BEGIN -check for near duplicate points - causes singularity in GMS
                                        d=GM.DistanceV(xs,ys,x[kx],y[ky])
                                        if np.min(d)<Dmin:
                                            if TextOutput:
                                                print("near duplicate point at distance: "+str(np.min(d))+" for node "+str(n))
                                            IncludePoint=False #END -check for near duplicate points
                                    zd=float(z.data)
                                    if zd < Zmin:#double check mask fail and bad fill value
                                            IncludePoint=False
    #                                if zd > Zmax:#exclude land values(should be optional)
    #                                        IncludePoint=False
                                    if IncludePoint: # first point
                                        zd=min(zd,Zmax) # trunkate interpolant to Zmax<=0
                                        xs.append(x[kx])
                                        ys.append(y[ky])
                                        zs.append(zd)
        else: # global scattered data, filetype[j] == 2
            fl=flnms[j]
            j=np.array( np.where(  np.array(np.abs(xp-x0) < SearchWidth[j] ) & np.array(np.abs(yp-y0) < SearchWidth[j] ) ) )
            j=list(itertools.chain.from_iterable(j))
            j=np.array(j).flatten()
            if j.size>0:
                j=np.array(j).flatten()
                x=x0[j]
                y=y0[j]
                zd=z0[j]
                jp=np.where( zd <= 0.)
                x=x[jp]
                y=y[jp]
                zd=zd[jp]
                if TextOutput:
                    print("number of points from "+fl+" is " +str(len(zd)) )
                xs.extend(x)
                ys.extend(y)
                zs.extend(zd)
    if TextOutput:      
        if n % 1 == 0:
            print("interpolating to:"+str(yp)+":"+str(xp))
            t1 = time.time()
            time_per_iter = (t1 - t0) / (n + 1)
            time_remaining = (nn - n) * time_per_iter / 60
            print("Progress: "+str(n+1)+" of "+ str(nn))
            print(f"Time remaining: {time_remaining:.2f} minutes")
            print(f"Average time per node: {time_per_iter:.2f} seconds")
            print("node: "+str(n)+" uses "+str(len(xs))+" data points")                        
    
    Npoints=len(xs)
    NumPoints[n]=Npoints

    #limit total number of points in interpolation (shouldn't come into effect much in practice)
    if Npoints > NpointsMax:
        if TextOutput:
            print("node "+str(n)+" has "+str(Npoints)+" data points. taking nearest: "+str(NpointsMax))
        xsp=np.zeros(NpointsMax)
        ysp=xsp
        zsp=xsp
        D=GM.Distance(xs,ys,xp,yp)
        for k in range(NpointsMax):
            ki=np.argmin(D)
            xsp[k]=xs[ki]
            ysp[k]=ys[ki]
            zsp[k]=zs[ki]
            D[ki]=float('inf')
        xs=xsp
        ys=ysp
        zs=zsp
        
##### Specify priori error statistics here#############################        
# Here some simple error statistics are given: These can be made dataset 
# dependent, depth dependent, etc.
# Assumptions regarding observation standard error:         
    ObsErr= max(1. , np.mean(np.abs(zs))/100. ) # 1 percent of local mean depth (m)
    VarObs=ObsErr**2 #(m^2)
#Specify assumptions regarding background variance:         
    VarBG=10.*VarObs #(m^2) Background variance = 10 observation error variance
##### Done specify priori error statistics ############################        

    
##### Different depth estimators ######################################        
    if TextOutput: 
        print("Number of obs = " +str(len(zs))  )
        print("std obs= " + str(np.sqrt(VarObs)) + " (m), std BG= " + str(np.sqrt(VarBG))+" (m)"  )

# Inverse distance estimator    
    ziID[n]=GM.InverseDistance(xs,ys,zs,xp,yp)

# Gauss-Markov smoothing with unkown mean (like ordinary kriging) 
    ziGMU[n], stdiGMU[n]=GM.GaussMarkovUnkMean(xs, ys, zs, xp, yp,LSp, VarObs,VarBG,True)
    if TextOutput: 
        print("GMU: est= "+str(ziGMU[n])+", err= "+str(stdiGMU[n]))

# Gauss-Markov smoothing with known mean (like simple kriging). Assumed mean is closest observation 
    ziGMN[n], stdiGMN[n]=GM.GaussMarkov(xs, ys, zs, xp, yp, LSp, VarObs,VarBG,"Nearest",True)
    if TextOutput: 
        print("GMN: est= "+str(ziGMN[n])+", err= "+str(stdiGMN[n]))

# Gauss-Markov smoothing with known mean (like simple kriging). Assumed mean is sample mean 
#   ziGMM[n]=GM.GaussMarkov(xs, ys, zs, xp, yp, LSp, VarObs,VarBG,"Mean",False)

# Gauss-Markov smoothing with known mean (like simple kriging). Assumed mean is 0 
#   ziGM0[n]=GM.GaussMarkov(xs, ys, zs, xp, yp, LSp, VarObs,VarBG,"Zero",False)

# Nearest neighbor estimator    
    D=GM.DistanceV(xs,ys,xp,yp)
    jc=np.argmin(D)
    ziClosest[n]=zs[jc]
    
###### Output estimates and other statistics to text files ###########        
#Inverse distance output
np.savetxt(floutID ,  ziID , fmt='%.6f', delimiter='\n')
#Gauss Markov smoothing with known mean output
np.savetxt(floutGMN,  ziGMN, fmt='%.6f', delimiter='\n')
#Std error (m) for Gauss Markov smoothing with known mean output
np.savetxt(floutGMNerr, stdiGMN, fmt='%.6f', delimiter='\n')
#Gauss Markov smoothing with unknown mean output
np.savetxt(floutGMU,  ziGMU, fmt='%.6f', delimiter='\n')
#Std error (m) for Gauss Markov smoothing with unknown mean output
np.savetxt(floutGMUerr,  stdiGMU, fmt='%.6f', delimiter='\n')
#Nearest neighbor estimate
np.savetxt(floutClosest,  ziClosest, fmt='%.6f', delimiter='\n')
#Number of bathymetry observations used for each node
np.savetxt(floutNpts, NumPoints, fmt='%i', delimiter='\n')
#Mesh length scale approcimation at nodes
np.savetxt(floutLLS, LocalLengthScale, fmt='%.6f', delimiter='\n')
    
    
