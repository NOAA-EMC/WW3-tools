function geom=pslg2geom(pslg)

geom.MSHID=3
geom.mshID= 'EUCLIDEAN-MESH'
geom.fileV=3
geom.point.coord=[pslg.x(:),pslg.y(:),0*pslg.x(:)];
[ne,two]=size(pslg.edges)
geom.edge2.index=[pslg.edges,zeros(ne,1)];
