function geom=pslg2geom(pslg)
%function geom=pslg2geom(pslg)
% convert Planar Straight Line Graph (pslg) structure with fields
%
% pslg.x : (nn x 1) x coordinates of nodes
% pslg.y : (nn x 1 )y coordinates of nodes
% pslg.edges : (nedges x 2) list of edges between nodes 
%
% to jigsaw geometry structure containing the same information
%

geom.MSHID=3
geom.mshID= 'EUCLIDEAN-MESH'
geom.fileV=3
geom.point.coord=[pslg.x(:),pslg.y(:),0*pslg.x(:)];
[ne,two]=size(pslg.edges)
geom.edge2.index=[pslg.edges,zeros(ne,1)];
