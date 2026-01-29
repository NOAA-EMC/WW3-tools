
import numpy as np
import netCDF4 as nc
import jigsawpy
import geopandas as gpd

def ReadWriteCoast(fl):

    gdf = gpd.read_file(fl)
    shp = gdf.shape
    n=shp[0]
    N=np.zeros((n,), dtype=int)
    xss=np.array([], dtype=np.single)
    yss=np.array([], dtype=np.single)

    for k in range(0, n-1):
        if np.mod(k,100)==0:
                print(str(k)+" "+str(n))
        seg=gdf.geometry[k]
        if seg.geom_type=='LineString':
            x,y = seg.coords.xy
        if seg.geom_type=='Polygon':
            bnd=seg.boundary
            x,y = bnd.coords.xy
        xp= np.array(x)
        yp= np.array(y)
        xss=np.append(xss,xp[:])
        yss=np.append(yss,yp[:])
      
        #x,y = seg.coords.xy
        #xp=numpy_array = np.array(x)
        #yp=numpy_array = np.array(y)
        ymp=np.mean(yp)
        dx=xp[:-1:] - xp[1::]
        dy=yp[:-1:] - yp[1::]
        lat2m=110574.
        lon2m=111320.*np.cos(ymp*np.pi/180.)
        d=np.sqrt( (dx*lon2m)**2 + (dy*lat2m)**2 )
        d=np.append(0,np.cumsum(d))
        nd=d.size-1
        npoints=int(np.ceil(d[nd]/dxS))
        di=np.linspace(d[0],d[nd],npoints)
        ni=di.size
        xi=np.zeros((ni,), dtype=np.single)
        yi=np.zeros((ni,), dtype=np.single)
        for j in range(0,ni):
            ks=np.where(  (d-di[j])**2  <  dxS**2  )
            xi[j]=np.mean(xp[ks])
            yi[j]=np.mean(yp[ks])
        xss=np.append(xss,xi[:])
        yss=np.append(yss,yi[:])
        if MakePlots:
            xo=np.append(xo,xp[:])
            yo=np.append(yo,yp[:])

    n=xss.size
    point=np.zeros((n,2), dtype=np.single)
    np.savetxt('CoastPoints.txt', (yss, xss), delimiter=' ')


    point[:,0] = xss[:]
    point[:,1] = yss[:]
    return point
