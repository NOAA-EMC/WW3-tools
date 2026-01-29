function g=CombineMesh(g0,g1);

nn=length(g0.x);
g.x=[g0.x(:);g1.x(:)];
g.y=[g0.y(:);g1.y(:)];
g.z=[g0.z(:);g1.z(:)];
g.e=[g0.e;nn+g1.e];

if isfield(g0,'bnd')
    if isfield(g1,'bnd')
        g.bnd=[g0.bnd;g1.bnd+nn];
    else
        g.bnd=g0.bnd
    end
else
    g.bnd=[];
end
    
