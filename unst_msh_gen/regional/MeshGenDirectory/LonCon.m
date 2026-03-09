function lon=LonCon(lon)
% function lon=LonCon(lon)
% The regional mesh generation assumes that the domain can be 
% represented on a plane (rather than the surface of a sphere). In 
% LonCon() we transform longitude to match longitude convention you
% want to use in your application. The transformation implemented for
% RWPS is used to avoid discontinuity at the international dateline 
% (-180,180) or Prime Meridian(360,0) etc.

% For most regional applications no transformation is nesesary so:

lon=lon; % Use for NWcoastal case and general regional mesh generation

% The transformation implemented for RWPS is used to avoid 
% discontinuity at the international dateline (-180,180) or Prime 
% Meridian(360,0) etc. The domain doesn't cross 90 deg E so we shift
% the discontinuity there.
% For RWPS specific mesh generation uncomment the line below:

%j=find(lon<90);lon(j)=lon(j)+360; % RWPS
