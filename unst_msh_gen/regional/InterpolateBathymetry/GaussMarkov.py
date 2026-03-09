
#
# Functions for Gauss Markov smoothing of bathymetry data.
# 
# Common input variables:
#   x,y,z are length nx vectors of observed depth (z) at point (x,y)
#   xi,yi is the points for z to be intepolated to.
#
# InverseDistance : weighted inverse distance (squared) estimator
#
# GaussMarkovUnkMean : Gauss Markov smoothing with unkown mean (like ordinary kriging)
#
# GaussMarkov :  Gauss Markov smoothing with kown mean (like simple kriging)
#   Choices for mean are : 
#       "Zero"      = 0 
#       "Mean"      = sample mean of z
#       "Median"    = sample median of z
#       "Nearest"   = value of z for closest point to (xi, yi)
#
# Two covariance functions are implemented:
#   covmod="Exp" for exponential decay covariance as a function of distance or
#   covmod="Sph" for spherical decay (finite suport)
#
# LengthScale is the (specified/ known) length scale in km at (xi, yi). 
# LengthScale is intended to be proportional to mesh lengthscake at (xi, yi).
#

import numpy as np

#covmod="Sph"
covmod="Exp"
CovExp=2 # Exponent in covariance funtion, used if covmod="Exp", 
IDExp=2 # exponent used in inverse distance estimator

deg2km=111.132954
deg2rad=np.pi/180

def InverseDistance(x,y,z,xi,yi):
    d=DistanceV(x,y,xi,yi)
    w=1.0/(d**IDExp)
    zi=np.dot( z, w )/np.sum( w )
    return zi

def GaussMarkov(x,y,z,xi,yi,LengthScale,Vobs,V,MeanTyp,CompErr):
    nx=len(x)
    f = np.zeros(nx)
    C = np.zeros((nx,nx))
    for j in range(nx):
        d = DistanceV(x,y,x[j],y[j])
        C[j,range(nx)]=CovarianceDistance(d,covmod,LengthScale,V)
        C[j,j]=C[j,j] + Vobs
    d = DistanceV(x,y,xi,yi)
    f = CovarianceDistance(d,covmod,LengthScale,V)
    detC = np.linalg.det(C)
    if np.isclose(detC, 0):
        zi=float("inf")
        print("Cov matrix is nearly singular. Check for double input points, etc.")
        if not CompErr:
            return zi
        else:
            erri=float("inf")
            return zi, erri
    else:
        match MeanTyp:
            case "Zero": # simple kriging type
                mu=0.
            case "Mean":
                mu=np.mean(z) # regional mean
            case "Median":
                mu=np.median(z) # regional median
            case "Nearest": # nearest value
                mu=z[np.argmax(f)] # max correlate(monotone)<==>closest point for process mean
            case _:
                mu=0.
        w = np.linalg.solve(C, f)
        zi=mu+np.dot(w,z-np.array(mu))
        if not CompErr:
            return zi
        else:
            erri=np.sqrt( V  -  np.dot(f,w)) 
            return zi, erri

def GaussMarkovUnkMean(x,y,z,xi,yi,LengthScale,Vobs,V,CompErr):
    nx=len(x)
    f = np.zeros(nx+1)
    C = np.zeros((nx+1,nx+1))
    for j in range(nx):
        d = DistanceV(x,y,x[j],y[j])
        C[j,range(nx)]=CovarianceDistance(d,covmod,LengthScale,V)
        C[j,j]=C[j,j] + Vobs
        C[nx,j]=1.
        C[j,nx]=1.
    d = DistanceV(x,y,xi,yi)
    f[range(nx)] = CovarianceDistance(d,covmod,LengthScale,V)
    C[nx,nx]=0
    f[nx]=1
    detC = np.linalg.det(C)
    if np.isclose(detC, 0):
        zi=float("inf")
        w=zi+f 
        print("Cov matrix is likely singular.")
    else:
        w = np.linalg.solve(C, f)
        zi= np.dot(w[range(nx)],z)
    if not CompErr:
        return zi
    else:
        erri=np.sqrt( V + np.dot( w[0:nx], np.matmul(C[0:nx,0:nx],w[0:nx]))  -2.* np.dot(w[0:nx],f[0:nx])   ) 
        return zi, erri

def DistanceV(x,y,x0,y0): #distance between list x,y and point x0,y0
    n=len(x)
    d=np.zeros(n)
    d=np.sqrt( 
        ( deg2km*np.cos( deg2rad*y0 )*(x-x0) )**2 + 
        ( deg2km*(y-y0))**2  ) 
    return d

#from geopy import distance
def Distance(x0,y0,x1,y1):
#    d = distance.great_circle(y0,x0,y1,x1).km
    d=np.sqrt( 
        ( deg2km*np.cos(deg2rad*(y0+y1)/2 )*(x1-x0) )**2 + 
        ( deg2km*(y1-y0))**2  ) 
    return d

def CovarianceDistance(d,model,ls,v):
    if model == "Exp":
        c = v*np.exp(- (d/ls )**CovExp )
    if model == "Sph":
        c = v*( 1. - 1.5*(d/ls) + .5*((d/ls)**3)   )
        c[np.where(d>ls)]=0.
    return c
        
