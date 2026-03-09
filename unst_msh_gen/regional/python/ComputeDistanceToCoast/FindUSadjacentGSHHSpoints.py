import numpy as np
import netCDF4 as nc
import jigsawpy
import geopandas as gpd
import time

from SmoothAndSubsampleCoastline import *

#compute points in the higher resolution GSHHS global coastline that are near the us coastline,
#defined in tigerline

#https://catalog.data.gov/dataset/tiger-line-shapefile-2019-nation-u-s-coastline-national-shapefile
#https://www.soest.hawaii.edu/pwessel/gshhg/

flUS="../RWPS/Data/us_coastline/tl_2023_us_coastline.shp"
GlobalCoastLineFile='../RWPS/Data/GlobalCoast/GSHHS_shp/f/GSHHS_f_L1.shp'

dxS=5000.       #(positive real)Smooth coastline to dxS meters length scale 
dxI=2500.       #(positive real)Interpolate smoothed coastline to dxI meters
pointsUS=SmoothAndSubsampleCoastline(flUS,dxS,dxI)

lonUS=pointsUS[:,0]
latUS=pointsUS[:,1]
#gdf = gpd.read_file(fl)
j=np.where( np.abs(lonUS+latUS)>= 0 ) # remove NaN's
lonUS=lonUS[j]
latUS=latUS[j]

pointsGlobal=SmoothAndSubsampleCoastline(GlobalCoastLineFile,dxS,dxI)
lonGlobal=pointsGlobal[:,0]
latGlobal=pointsGlobal[:,1]

np.savetxt('GlobalCoastPoints.txt', (latGlobal, lonGlobal), delimiter=' ')

j=np.where( np.abs(lonUS+latUS)>= 0 ) # remove NaN's
lonUS=lonUS[j]
latUS=latUS[j]

j=np.where( np.abs(lonGlobal+latGlobal)>= 0 ) # remove NaN's
lonGlobal=lonGlobal[j]
latGlobal=latGlobal[j]

npGlobal=lonGlobal.size


t0=time.time()
D=np.zeros((npGlobal),dtype=np.single)
lat2m=np.single(110574.)
i=complex(0,1)
for k  in range(0,npGlobal-1):
    lon2m=np.single(111320.*np.cos(latGlobal[k]*np.pi/180.))
    DLON=np.mod( lonGlobal[k] - lonUS , 360. )
    DLAT=        latGlobal[k] - latUS
    ldE=np.min (  np.abs(           DLON*lon2m + i*DLAT*lat2m  ) )#Distance to east
    ldW=np.min (  np.abs(  (360. - DLON)*lon2m + i*DLAT*lat2m  ) )#Distance to west
    D[k]=min(ldE,ldW) # min of east, west distances
    t1=time.time()
    tps=(t1-t0)/max(k,1)
    tre=tps*(npGlobal-k)
    if np.mod(k,10000)==0:
        print( str(k) + " of " + str(npGlobal)+" percent:"+str(100*k/npGlobal) )
#        print("estimated hours remaining: "+ str( tre /3600. ) )
        print("estimated minutes remaining: "+ str( tre /60. ) )

np.savetxt('USCoastPoints.txt', (latUS, lonUS), delimiter=' ')

thrsh=111000.
j=np.where( D <= thrsh ) # remove NaN's
np.savetxt('USCoastPointsFromGSHHS.txt', (latGlobal[j], lonGlobal[j]), delimiter=' ')

lonUSX=np.append(lonUS,lonGlobal[j])
latUSX=np.append(latUS,latGlobal[j])
np.savetxt('USCoastPointsWithGSHHS.txt', (latUSX, lonUSX), delimiter=' ')


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#add points to off shore banks we want to refine here we have Georges Bank(GB)
#and banks around the Bahamas(FB) as well as Penguin bank of Molakai, HI.
#This artificially increases resolution for these features
lonFB=[  -78.6485,  -78.6999,  -78.3401,  -76.4383,  -79.9849]
latFB=[   27.0490,   25.4556,   23.6567,   22.9114,   23.7338]
#lonFB= [ -79.9909,  -78.1963,  -78.3957]
#latFB=  [23.7978,   26.8928,   24.1530]
lonUSX=np.append(lonUSX,lonFB)
latUSX=np.append(latUSX,latFB)
lonGB=[  -67.4977,  -68.0673]
latGB=[   41.6880,   41.1443]
#lonGB = -67.4517
#latGB =   41.3567
lonUSX=np.append(lonUSX,lonGB)
latUSX=np.append(latUSX,latGB)
#penguin bank HI 50m deep- probably not important in waves
lonPB =[ -157.6623, -157.5066]
latPB=[   20.9413,   21.0486]
lonUSX=np.append(lonUSX,lonPB)
latUSX=np.append(latUSX,latPB)

np.savetxt('USCoastPointsWithGSHHSandBanks.txt', (latUSX, lonUSX), delimiter=' ')

