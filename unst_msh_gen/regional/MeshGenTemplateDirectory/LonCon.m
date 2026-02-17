function lon=LonCon(lon)
% function lon=LonCon(lon)
% Transform longitude to match longitude convention you want to use
% here.  This is used to avoid discontinuity at the international;
% dateline etc.

 
lon=lon;

%For RWPS meshes spanning 180W,180E but not 90E
% use below:
%j=find(lon<90);lon(j)=lon(j)+360; % RWPS

