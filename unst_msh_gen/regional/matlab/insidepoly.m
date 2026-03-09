function [inpoly,onpoly0]=insidepoly( xi,yi,xpoly,ypoly);
inpoly=inpoly2([xi(:),yi(:)],[xpoly(:),ypoly(:)]);

%Dummy output for onpoly
if nargout>1,
    onpoly0=0*inpoly;
end
%INPOLY2 compute "points-in-polygon" queries.  
%   [STAT] = INPOLY2(VERT,NODE,EDGE) returns the "inside/ou-
%   tside" status for a set of vertices VERT and a polygon 
%   {NODE,EDGE} embedded in a two-dimensional plane. General
%   non-convex and multiply-connected polygonal regions can 
%   be handled. VERT is an N-by-2 array of XY coordinates to 
%   be tested. STAT is an associated N-by-1 logical array,
%   with STAT(II) = TRUE if VERT(II,:) is an interior point.
