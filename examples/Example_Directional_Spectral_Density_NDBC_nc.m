clear all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is an example convenrt ndbc spectral file into%
% WW3 netcdf format                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Ali Abdolali Feb 2023 ali.abdolali@noaa.gov          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add path
addpath ../matlab_tools
%% input data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfreq=35;                  % number of frequencies
nDir=36;                   % number of Directions
Dir0=5;                    % first direction (deg)
inc=1.1;                   % frequency increment
f0=0.038; 
filename='bnd_ndbc_ww3.nc';% name of netcdf file (boundary)
testcase='DUCK';           % test vase
pointID='44100';           % id (length should be 16)
Lat=34.71;                 % latitude of BC (deg)
Lon=-72.248;               % longitude of BC (deg)
dpt=26;                    % Depth of BC (m)
wndspd=0;                  % wind speed (m/s)
wnddir=270;                % wind direction (deg)
curspd=0;                  % current velocity (m/s)
curdir=270;                % current direction (deg)
coordinate='spherical';    % coordinate : Spherical, cartesian
ncf='44100w9999.nc';       % input netcdf file from https://dods.ndbc.noaa.gov/
visualize='true';          % options:true, false- this requires polarPcolor function and can be obtained from ...
                           % https://www.mathworks.com/matlabcentral/fileexchange/49040-pcolor-in-polar-coordinates
deltatheta=10;             % DeltaDir
theta0=0;                  % first Dir
%----------------------------------------------------------%
% if visualize='true', user can define the time frame to visualize
t1=(datenum('20221201 000000','yyyymmdd HHMMSS')); % first timestep 
t2=(datenum('20221202 000000','yyyymmdd HHMMSS')); % last timestep 
dt=1/24;                      %delta t (day)
%----------------------------------------------------------%
% frequency (log)
Omega(1)=2*pi*f0;
for i=2:nfreq
Omega(i)=inc*Omega(i-1);
end
freq=Omega/2/pi;
%----------------------------------------------------------%
display (['Reading ', ncf,' ...'])
% read input netcdf file (directional spectral density file)
[SWDEN] = swden_ndbc_read(ncf,deltatheta,theta0,freq);
%----------------------------------------------------------%
% The boundary requres coordinates (lat, lon) depth, wind speed
% wind direction, current speed, current direction for each time step.
% Here wind and current are 0.
pointID = [pointID,repmat(' ', [1, 16-strlength(pointID)])];
Lon=Lon*ones(1,length(SWDEN.Int.time));
Lat=Lat*ones(1,length(SWDEN.Int.time));
curspd=curspd*ones(1,length(SWDEN.Int.time));
curdir=curdir*ones(1,length(SWDEN.Int.time));
wndspd=wndspd*ones(1,length(SWDEN.Int.time));
wnddir=wnddir*ones(1,length(SWDEN.Int.time));
dpt=dpt*ones(1,length(SWDEN.Int.time));
time(1,:)=SWDEN.Int.time;
dir=pi*SWDEN.Int.Dir/180; %radian
dir0=SWDEN.Int.Dir;%degree
EFTH(:,:,1,:)=SWDEN.Int.DENS; %directional spectral density time series
%----------------------------------------------------------%
display (['Generating ', filename,' ...'])
%dump into netcdf
[filename] = write_directional_spectra_nc(filename,testcase,...
            pointID,Lat,Lon,dpt,wndspd,wnddir,curspd,curdir,time,...
            SWDEN.Int.f,dir0,EFTH,coordinate);
        
 
%%
%----------------------------------------------------------%        
tf=strcmp(visualize,'true');                 
if tf==1
display (['Plot Generation ', [filename,'.gif'],' ...'])
width=1200;  % Width of figure for movie [pixels]
height=500;  % Height of figure of movie [pixels]
left=200;     % Left margin between figure and screen edge [pixels]
bottom=200;  % Bottom margin between figure and screen edge [pixels]

DENSOrig(:,:,:)=SWDEN.Orig.DENS(:,:,:);
DENSInt(:,:,:)=SWDEN.Int.DENS(:,:,:);
t=t1:dt:t2;
h=figure;
set(gcf,'Position', [left bottom width height])

for k=1:1:length(t)
    clf
 subplot(1,2,1)
 
 [it]=find(abs(SWDEN.Orig.time-t(k))==min((abs(SWDEN.Orig.time-t(k)))));
 clear DENSS
 DENSS(:,:)=DENSOrig(:,:,it(1),:);
 DENSStmp=[DENSS;DENSS(1,:)];
 [~,cc] = polarPcolor(SWDEN.Orig.f,[SWDEN.Orig.Dir SWDEN.Orig.Dir(1)],...
     DENSStmp,'Nspokes',36,'ncolor',10,'labelR','f (Hz)','Rscale','log');
 ylabel(cc,['Obs (Original); Date = ',datestr(SWDEN.Orig.time(it(1)))],'FontSize',14);
caxis([0 nanmax(DENSS(:))])
 subplot(1,2,2)
 clear DENSS
 DENSS(:,:)=DENSInt(:,:,it(1),:);
 DENSStmp=[DENSS;DENSS(1,:)];
 [~,cc] = polarPcolor(SWDEN.Int.f,[SWDEN.Int.Dir SWDEN.Int.Dir(1)],...
     DENSStmp,'Nspokes',36,'ncolor',10,'labelR','f (Hz)','Rscale','log');
 ylabel(cc,['Obs (Interpolated); Date = ',datestr(SWDEN.Orig.time(it(1)))],'FontSize',14);
caxis([0 nanmax(DENSS(:))])

  frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    if k == 1
        imwrite(imind,cm,[filename,'.gif'],'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,[filename,'.gif'],'gif','WriteMode','append');
    end   
end 
end
%----------------------------------------------------------%         
display (['Finished'])
