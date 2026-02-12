import os
import numpy as np
import netCDF4 as nc
import jigsawpy
import geopandas as gpd
import time
from SmoothAndSubsampleCoastline import *

MakePlots=False

if MakePlots:
    import plotly.express as px
    import plotly.graph_objects as go

#(positive integer)Subsample bathymetry coordinates to Nsubsamble
Nsubsample=2 #86 hours

Nsubsample=5 #3.75 hours

flout="DistToUSCoast.msh"
CoastInput=1
if CoastInput==0:
    fl="../RWPS/Data/us_coastline/tl_2023_us_coastline.shp"
    #gdf = gpd.read_file("../RWPS/Data/us_coastline/tl_2023_us_coastline.shp")
    dxS=5000.       #(positive real)Smooth coastline to dxS meters length scale 
    dxI=2500.       #(positive real)Interpolate smoothed coastline to dxI meters
    points=SmoothAndSubsampleCoastline(fl,dxS,dxI)
    lonUS=points[:,0]
    latUS=points[:,1]

if CoastInput==1: # Use precomputed smoothed & subsampled coast points from FindUSadjacentGSHHSpoints.py
    fl="./USCoastPointsWithGSHHSandBanks.txt"
    points=np.loadtxt(fl, delimiter=' ')
    lonUS=points[1,:]
    latUS=points[0,:]

#gdf = gpd.read_file(fl)
j=np.where( np.abs(lonUS+latUS)>= 0 ) # remove NaN's
lonUS=lonUS[j]
latUS=latUS[j]

np.savetxt('CoastPoints.txt', (lonUS, latUS), delimiter=' ')
    
print("number of coastline points = " + str(lonUS.size))
print("This number has a strong effect on run time")
    
# load topo file to build distance function on the grid of
data = nc.Dataset("../RWPS/Data/RTopo_2_0_4_GEBCO_v2023_60sec_pixel.nc","r")
xlon = np.asarray(data["lon"][:])
ylat = np.asarray(data["lat"][:])
#elev = np.asarray(data["bed_elevation"][:]) + np.asarray(data["ice_thickness"][:])

if MakePlots:
    fig = px.scatter(x=xlon, y=ylat)
    fig.show(renderer='browser')

xmid = 0.5 * (xlon[:-1:] + xlon[1::])
ymid = 0.5 * (ylat[:-1:] + ylat[1::])

xmid=xmid[0 :: Nsubsample]
ymid=ymid[0 :: Nsubsample]

nx=xmid.size
ny=ymid.size
npoints=lonUS.size

D=np.zeros((ny,nx),dtype=np.single)
lat2m=np.single(110574.)
i=complex(0,1)

t0=time.time()
D=np.zeros((ny,nx),dtype=np.single)
lat2m=np.single(110574.)
i=complex(0,1)
DLON=np.mod( np.array([xmid] * npoints) - np.array([lonUS] * nx).T, 360. )
for k  in range(0,ny-1):
    print( str(k) + " of " + str(ny)+" percent:"+str(100*k/ny) )
    lon2m=np.single(111320.*np.cos(ylat[k]*np.pi/180.))
    DLAT=(ymid[k]-latUS)
    DLATnx=np.array([DLAT] * nx).T
    ldE=np.min (  np.abs(           DLON*lon2m + i*DLATnx*lat2m  ), axis=0 )#Distance to east
    ldW=np.min (  np.abs(  (360. - DLON)*lon2m + i*DLATnx*lat2m  ), axis=0 )#Distance to west
    ld=np.array([ldE,ldW])
    D[k,:]=np.min(ld,axis=0) # min of east, west distances
    t1=time.time()
    tps=(t1-t0)/max(k,1)
    tre=tps*(ny-k)
    print("estimated hours remaining: "+ str( tre /3600. ) )

#output to jigsaw .msh format
dist = jigsawpy.jigsaw_msh_t() 
dist.mshID = "ellipsoid-grid"
dist.radii = np.full( +3, +6371.0, dtype=jigsawpy.jigsaw_msh_t.REALS_t)
dist.xgrid = xmid
dist.ygrid = ymid 
dist.value = D
jigsawpy.savemsh(flout, dist)


#make some simple plots for sanity check
if MakePlots:
    fig = px.imshow(D)
    fig.show(renderer='browser')

    fig1 = go.Figure(data=go.Heatmap(z=D/1000,x=xmid,y=ymid))

    fig1.update_layout(
        title=dict(text='Distance to US coastline (km)'),
        xaxis_nticks=36)
    fig1.show(renderer='browser')
    fig1.add_trace(
        go.Scatter(
            x=lonUS,
            y=latUS,
            mode='markers',
            showlegend=False)
        )
    fig1.show(renderer='browser')

if MakePlots:

    fig2 = px.imshow(D/1000,labels=dict(x="longitude", y="latitude", color="distance"),
                    x=xmid,
                    y=ymid
                    )

    fig2.update_layout( title=dict(text='Distance to US coastline (km)') )
    fig2.show(renderer='browser')
    fig2.show(renderer='browser')
    fig2.add_trace(
        go.Scatter(
            x=lonUS,
            y=latUS,
            mode='markers',
            showlegend=False)
        )
    fig2.show(renderer='browser')


    fig = px.scatter(x=lonUS, y=latUS)
    fig.add_trace(px.scatter(x=lonUS[1], y=latUS[1], mode='markers', marker=dict(color='red'), name='Second Trace'))
    #fig.add_trace(go.Scatter(x=x2, y=y2, mode='markers', marker=dict(color='red'), name='Second Trace'))
    fig.show(renderer='browser')
