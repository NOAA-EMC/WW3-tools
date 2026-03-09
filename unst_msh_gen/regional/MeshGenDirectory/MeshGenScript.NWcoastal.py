
import os
import argparse
import time
import numpy as np
import netCDF4 as nc

import jigsawpy

from scipy.interpolate import RegularGridInterpolator

#-------------------Input Files----------------------------------------

#GSSH coastline 
#jigsaw .msh format Planer Straight Line Graph defining mesh outer boundary and coastline
PSLGFile="NWcoastal.PSLG.msh"
#jigsaw gridded .msh format Distance to taget poings
DistanceToCoastFile="DFun.NWcoastal.PSLG.msh"
#jigsaw gridded .msh format topography on same grid as distance
TopographyFile="Topo.DFun.NWcoastal.PSLG.msh"

# directory to write output files to
OutDir='NWcoastal/'

#-------------------Paramter Inputs------------------------------------
#parameters for specifying resolution
d0=15000.
d1=300000.
beta=2000.
Smin=0.25
Smax=7.5

#-------------------Main Program---------------------------------------

# Create the output directory------------------------------------------
try:
    os.mkdir(OutDir)
    print(f"Directory '{OutDir}' created successfully.")
except FileExistsError:
    print(f"Directory '{OutDir}' already exists. Proceeding ...")
except PermissionError:
    print(f"Permission denied: Unable to create '{OutDir}'.")

# Setup jigsaw structures----------------------------------------------
opts = jigsawpy.jigsaw_jig_t()
topo = jigsawpy.jigsaw_msh_t() #topographic database
dist = jigsawpy.jigsaw_msh_t() #distance to US coastline in meters
geom = jigsawpy.jigsaw_msh_t()
mesh = jigsawpy.jigsaw_msh_t()
meshR3 = jigsawpy.jigsaw_msh_t()
hmat = jigsawpy.jigsaw_msh_t()
proj = jigsawpy.jigsaw_prj_t()

opts.geom_file = "geom.msh"
opts.jcfg_file = "jcfg.jig"
opts.mesh_file = "mesh.msh"
opts.hfun_file = "spac.msh"

# load input data------------------------------------------------------
jigsawpy.loadmsh(PSLGFile, geom)
jigsawpy.loadmsh(DistanceToCoastFile, dist)
jigsawpy.loadmsh(TopographyFile, topo)

# truncate data to bounding rectangle of input PSLG--------------------
xmin = np.min( geom.point["coord"][:, 0])
ymin = np.min( geom.point["coord"][:, 1])
xmax = np.max( geom.point["coord"][:, 0])
ymax = np.max( geom.point["coord"][:, 1])

tv = topo.value
dv = dist.value

xmsk = np.logical_and( topo.xgrid > xmin , topo.xgrid < xmax )
ymsk = np.logical_and( topo.ygrid > ymin , topo.ygrid < ymax )

tv = tv[:, xmsk]
tv = tv[ymsk, :]
dv = dv[:, xmsk]
dv = dv[ymsk, :]

# define spatial resolution shape function----------------------------- 
W=np.exp(-  (np.abs(dv-d0) / d1) - (np.abs(tv) / beta)  )
W[ np.where(dv < d0 ) ] = 1.

# build resolution specification gridded data structure, hmat----------
hmat=topo
hmat.value=Smax - (Smax - Smin)*W

hmat.mshID = "ellipsoid-grid"
hmat.radii = np.full( +3, +6371.0, dtype=jigsawpy.jigsaw_msh_t.REALS_t)

hmat.xgrid = hmat.xgrid[xmsk] * np.pi / 180.
hmat.ygrid = hmat.ygrid[ymsk] * np.pi / 180.

hmat.value = np.maximum(hmat.value, Smin)#should not be nescesary
hmat.value = np.minimum(hmat.value, Smax)

hmat.slope = np.full( hmat.value.shape, +0.1500 , dtype=jigsawpy.jigsaw_msh_t.REALS_t)

#------------------------------------ do stereographic proj.
geom.point["coord"][:, :] *= np.pi / 180.

proj.prjID = 'stereographic'
proj.radii = +6.371E+003
proj.xbase = +0.500 * (xmin + xmax) * np.pi / 180.
proj.ybase = +0.500 * (ymin + ymax) * np.pi / 180.

jigsawpy.savemsh(OutDir+"HFUN.msh", hmat)

jigsawpy.project(geom, proj, "fwd")
jigsawpy.project(hmat, proj, "fwd")

jigsawpy.savemsh(opts.geom_file, geom)
jigsawpy.savemsh(opts.hfun_file, hmat)

# save hmat------------------------------------------------------------
jigsawpy.savemsh(OutDir+"HFUNproj0.msh", hmat)

#smooth hmat
jigsawpy.cmd.marche(opts, hmat)

# save smoothed hmat---------------------------------------------------
jigsawpy.savemsh(OutDir+"HFUNproj1.msh", hmat)

# make mesh using JIGSAW-----------------------------------------------
opts.hfun_scal = "absolute"
opts.hfun_hmax = float("inf")       # null HFUN limits
opts.hfun_hmin = float(+0.00)
#opts.hfun_hmax = 1.25*Smax       # Unintended effects, better off null 
#opts.hfun_hmin = .5*Smin

opts.mesh_dims = +2                 # 2-dim. simplexes
opts.mesh_eps1 = +1.

#opts.mesh_top1 = "true" !!!No convergece


ttic = time.time()

jigsawpy.cmd.jigsaw(opts, mesh)

ttoc = time.time()

print("CPUSEC =", (ttoc - ttic))

# save mesh in JIGSAW native projection based on input PSLG------------
jigsawpy.savemsh(OutDir+"/RWPS.PROJ.msh",mesh)

# compute costa functions and save
cost = jigsawpy.triscr2(mesh.point["coord"],mesh.tria3["index"])
np.savetxt(OutDir+"TriScr2.txt",cost,"%f")
print("TRISCR =", np.min(cost), np.mean(cost))

cost = jigsawpy.pwrscr2(mesh.point["coord"], mesh.power, mesh.tria3["index"])
np.savetxt(OutDir+"PwrScr2.txt",cost,"%f")
print("PWRSCR =", np.min(cost), np.mean(cost))

tbad = jigsawpy.centre2(mesh.point["coord"],mesh.power,mesh.tria3["index"])
print("OBTUSE =",+np.count_nonzero(np.logical_not(tbad)))


# project mesh nodes to radian lat, lon
jigsawpy.project(mesh, proj, "inv") # This used to work
#jigsawpy.savemsh("RWPS.radian.msh",mesh)

# transform mesh nodes to degree lat, lon
mesh.point["coord"][:, :] = mesh.point["coord"][:, :]*180. / np.pi

# save mesh in degree lat, lon system
jigsawpy.savemsh(OutDir+"RWPS.LL.msh",mesh)

# create jigsaw R3 mesh on global surface and save to 
S2=mesh.point["coord"][:,[0,1]]
S2=S2*np.pi/180.

R3=jigsawpy.S2toR3(mesh.radii,S2)
#np.savetxt("R3.txt",R3," %f ") # save 3D nodes

meshR3 = jigsawpy.jigsaw_msh_t()
mesh.mshID = 'ellipsoid-mesh'
meshR3.tria3=mesh.tria3
meshR3.ndims=3
#make 3D coordinates 
nd=R3.shape
meshR3.vert3 = np.zeros(nd[0], dtype=mesh.VERT3_t)
meshR3.vert3["coord"] = R3
jigsawpy.savemsh(OutDir+"RWPS.R3.msh",meshR3)


# Now apply filters and output in WW3 form
# import FilterRoutines.py

from FilterRoutinesNM import *

#replace mesh with R3 mesh, mesh->mesh R3
jigsawpy.loadmsh(OutDir+"RWPS.R3.msh", mesh) # uncomment if starting here 

opts.geom_file = "geom.msh"  #saves the geometry info for jigsaw
opts.jcfg_file = "opts.jig"  #jigsaw ctlr file
   
    
geom.mshID = "ellipsoid-mesh"
geom.radii = np.full(3, 6.371E+003, dtype=geom.REALS_t)

jigsawpy.savemsh(opts.geom_file, geom)

inject_dem()
filter_ocn()
    
jigsawpy.savemsh(OutDir+"RWPS.F.R3.msh", mesh)

# viz. in eg. paraview
jigsawpy.savevtk(OutDir+"test.vtk", mesh)

# convert to lon lat    
point = mesh.point["coord"]
point = jigsawpy.R3toS2(geom.radii, point) 
point*= 180. / np.pi
 
depth = np.reshape(-1*mesh.value, (mesh.value.size, 1))
depth[depth <= 0] = 2
point = np.hstack((point, depth))  # append elev. as 3rd coord.
cells = [("triangle", mesh.tria3["index"])]
tri_data=cells[0][1]+1

#put coordinates in non standard format to avoid international date line
lon=point[:,0]
lon[np.where(lon>90)]=lon[np.where(lon>90)]-360
point[:,0]=lon

write_gmsh_mesh(OutDir+"RWPS.ww3", point, tri_data)

mesh.point["coord"]=point
jigsawpy.savemsh(OutDir+"RWPS.F.LLH.msh", mesh)


#write final mesh in jigsaw .msh format
meshR2 = jigsawpy.jigsaw_msh_t()
#make 2D coordinates 
nd=point.shape
meshR2.ndims=2
meshR2.vert2 = np.zeros(nd[0], dtype=mesh.VERT2_t)
meshR2.vert2["coord"] = point[:,[0,1]]
meshR2.tria3=mesh.tria3
meshR2.mshID=mesh.mshID

jigsawpy.savemsh(OutDir+"RWPS.F.LL.msh", meshR2)
    
