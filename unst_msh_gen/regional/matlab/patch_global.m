function ph=patch_global(lon,lat,f,e)
%function ph=patch_global(lon,lat,f,e)
% script to make a patch graphic object showing field f structure for mesh with coordinates lon,lat and element matrix e.
% elements crossing dateline are eliminated in the plot

    dx=max(lon(e))-min(lon(e));
    j=find(dx<270);
    ep=e(:,j);
    
    ph=patch(lon(ep),lat(ep),f(ep));
