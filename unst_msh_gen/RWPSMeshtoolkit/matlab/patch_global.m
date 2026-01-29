function ph=patch_global(lon,lat,f,e)

    dx=max(lon(e))-min(lon(e));
    j=find(dx<270);
    ep=e(:,j);
    
    ph=patch(lon(ep),lat(ep),f(ep));