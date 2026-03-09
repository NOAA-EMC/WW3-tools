
import numpy as np
import netCDF4 as nc
import jigsawpy
import geopandas as gpd

MakePlots=False
MakePlots=True

if MakePlots:
    import plotly.express as px
    import plotly.graph_objects as go


def SmoothAndSubsampleCoastlineP(fl,dxS,dxI):

    gdf = gpd.read_file(fl)
    shp = gdf.shape
    n=shp[0]
    N=np.zeros((n,), dtype=int)
    for k in range(1, n):
        if np.mod(k,1000)==0:
            print("Loading Coastline: "+str(k)+" of "+str(n) )
        try:
            seg=gdf.geometry[k-1]
            bnd=seg.boundary
            x,y = bnd.coords.xy
            N[k-1]=x.size
        except:
            print("Error reading segment # "+str(k)) 
            N[k-1]=0

    k = np.argmax(N)
    m = np.max(N)
    
    print("smoothing "+str(n)+" coast segments. longest segment at "+str(k)+" of length "+str(m))
    #gdfS=gdf
    xss=np.array([], dtype=np.single)
    yss=np.array([], dtype=np.single)
    
    if MakePlots:
        xo=np.array([], dtype=np.single)
        yo=np.array([], dtype=np.single)
    
    for k in range(0, n-1):
        if np.mod(k,100)==0:
                print(str(k)+" "+str(n))
        seg=gdf.geometry[k]
        bnd=seg.boundary
        x,y = bnd.coords.xy
        xp=numpy_array = np.array(x)
        yp=numpy_array = np.array(y)
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

    point[:,0] = xss[:]
    point[:,1] = yss[:]
    if MakePlots:
#        fig = px.scatter(x=xss, y=yss, mode='markers', marker=dict(color='blue'), name='smoothed + subsampled coastline')
#        fig = px.scatter(x=xss, y=yss)
#        fig.add_trace(go.Scatter(x=xo, y=yo, mode='markers', marker=dict(color='red'), name='input coastline'))
#        fig.show(renderer='browser')
        xo=xo[0 :: 100]
        yo=yo[0 :: 100]
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=xss, y=yss,mode='markers',name='smoothed'))#,color='blue'))
        # go.Scatter can't handle more than ~1 mm points-  so this may not work
        fig.add_trace(go.Scatter(x=xo , y=yo ,mode='markers',name='origonal'))#,color='red'))
        fig.show()
    
    return point
