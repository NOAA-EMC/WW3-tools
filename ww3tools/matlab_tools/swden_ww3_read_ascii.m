function [SWDEN] = swden_ww3_read_ascii(specfile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function reads WW3 directional spectral density file %
% .spec asill format                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Ali Abdolali Feb 2023 ali.abdolali@noaa.gov          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  input data %--------------------------------------------%
% specfile: name of .spec file
%  output data %--------------------------------------------%
% buoy_name: list of buoy
% time (Matlab time)
% longitude of points (Degree)
% latitude of points (Degree)
% f: frequency (Hz)
% Dir: direction (degree)
% DENS: Directional Spectral Density [freq,direction,stations,time]
% SPEC: Spectral Density [freq,stations,time]
% Hs: Significant Wave Heigth (m)
% Fp: Peak Freq (Hz)
% mwvdir: mean wave direction (degrees)
% cur: current speed (m/s)
% curdir: current direction (deg)
% wnd: wind speed (m/s)
% wnddir: wind direction (deg)
% dpt: depth (m)

%----------------------------------------------------------%
m=0;
%read and put all lines into cells
fid = fopen(specfile,'rt');
tline = fgetl(fid);
while ischar(tline)
m=m+1;
     file{m}=(tline);
    tline = fgetl(fid);
end
fclose(fid);
%----------------------------------------------------------%
freq=[];
Dir=[];
%'---------------------'    nfreq nDir   nstation '------------------------'
%'WAVEWATCH III SPECTRA'     50    36     1 'spectral resolution for points'
C = textscan(file{1},'%s %s %s %n %n %n %s %s %s %s');
nfreq=C{4};
nDir=C{5};
nStation=C{6};
  i=1;
    while length(freq)<nfreq
        i=i+1;
        A=cell2mat(textscan(file{i}, '%n', 'delimiter', '\n', 'whitespace', ' '));
    freq=[freq;A];
    end
%----------------------------------------------------------%    
  while length(Dir)<nDir
        i=i+1;
        A=cell2mat(textscan(file{i}, '%n', 'delimiter', '\n', 'whitespace', ' '));
    Dir=[Dir;A];
  end
%----------------------------------------------------------%
%go through time loop
  m=1;
   while i<length(file)
  i=i+1;
      time(m,1) = datenum(file{i},'yyyymmdd HHMMSS');
   i=i+1;
% buoyname  '  lat    long       dpt     wnd  wnddir cur  curdir   
%'46069     '  33.65-120.20     919.8   5.10  30.0   0.38 356.0
    C=textscan(file{i}, '%s %s %n %n %n %n %n %n %n %*[^\n]');
    buoy_name=C{1};
    lat(1,m)=C{3};
    lon(1,m)=C{4};
    dpt(1,m)=C{5};
    wnd(1,m)=C{6};
    wnddir(1,m)=C{7};
    cur(1,m)=C{8};
    curdir(1,m)=C{9};
%----------------------------------------------------------%    
  dens=[];
    while length(dens)<nDir*nfreq
        i=i+1;
        A=cell2mat(textscan(file{i}, '%n', 'delimiter', '\n', 'whitespace', ''));
    dens=[dens; A];
    end
  denstmp(1,:)=dens;
  DENS(:,:,nStation,m)=fliplr(reshape(dens,[nDir,nfreq]));
   m=m+1;
   end
  
   
   
SWDEN.buoy_name=cellstr(buoy_name);
SWDEN.time=time;
SWDEN.latitude=lat;
SWDEN.longitude=lon;
SWDEN.f=freq;
SWDEN.Dir=round(Dir*180/pi);%precision is not enough in spec
SWDEN.DENS=DENS;
SWDEN.cur=cur;
SWDEN.curdir=curdir;
SWDEN.wnd=wnd;
SWDEN.wnddir=wnddir;
SWDEN.dpt=dpt;
DeltaDir=abs(SWDEN.Dir(2)-SWDEN.Dir(1)); 
%DeltaDir
DENS=SWDEN.DENS;
DENS=[DENS; DENS(1,:,:,:)];
[DD,FF,KK,TT]=size(SWDEN.DENS);

DENSCOS=cosd([SWDEN.Dir;SWDEN.Dir(1)]).*DENS;
DENSSIN=sind([SWDEN.Dir;SWDEN.Dir(1)]).*DENS;

for k=1:KK
for ii=1:FF
SPEC(ii,k,:) = abs(trapz(0:pi*DeltaDir/180:2*pi,DENS(:,ii,k,:)));
SPECA(ii,k,:) = (trapz(0:pi*DeltaDir/180:2*pi,DENSCOS(:,ii,k,:)));
SPECB(ii,k,:) = (trapz(0:pi*DeltaDir/180:2*pi,DENSSIN(:,ii,k,:)));
end
end
SWDEN.SPEC=SPEC;
for k=1:KK
 Hs(k,:)=4.004*sqrt(abs(trapz(2*pi*SWDEN.f,SPEC(:,k,:))))/sqrt(2*pi);
 AA(k,:)=((trapz(2*pi*SWDEN.f,SPECA(:,k,:))));
 BB(k,:)=((trapz(2*pi*SWDEN.f,SPECB(:,k,:))));
end
 SWDEN.Hs=Hs;
 SWDEN.mwvdir=(180*atan2(BB,AA)/pi)+180;

for k=1:KK
    clear SPECmax
    clear SPECtmp
    SPECtmp(:,:)=SPEC(:,k,:);
    SPECmax=nanmax(SPECtmp);
for i=1:length(SWDEN.time)
    clear i1
    clear j1
    [i1,j1]=find(SPECtmp(:,i)==SPECmax(i));
    Fp(k,i)=nanmean(SWDEN.f(i1));
end
SWDEN.Fp=Fp;    
end
