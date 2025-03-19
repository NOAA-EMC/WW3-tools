clear;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is an example to generate Jonswap spectrum in %
% WW3 netcdf format                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Ali Abdolali Feb 2023 ali.abdolali@noaa.gov          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add path
addpath ../matlab_tools
%% input data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hs=1;                      % significant wave height (m)                  
Tp=13;                     % peak period (s)                              
MeanDir=100;                % Mean Wave direction
nfreq=35;                  % number of frequencies
nDir=36;                   % number of Directions
Dir0=5;                    % first direction (deg)
inc=1.1;                   % frequency increment
f0=0.038;                  % first freq (Hz)
spread=2;                  % directional spread (times DeltaDir)
filename='boundary.nc';    % name of netcdf file (boundary)
testcase='DUCK';           % test vase
pointID='wavemaker';       % id
Lat=34.71;                 % latitude of BC (deg)
Lon=-72.248;               % longitude of BC (deg)
dpt=26;                    % Depth of BC (m)
wndspd=0;                  % wind speed (m/s)
wnddir=270;                % wind direction (deg)
curspd=0;                  % current velocity (m/s)
curdir=270;                % current direction (deg)
coordinate='spherical';    % coordinate : Spherical, cartesian
visualize='true';          % options:true, false- this requires polarPcolor function and can be obtained from ...
                           % https://www.mathworks.com/matlabcentral/fileexchange/49040-pcolor-in-polar-coordinates
%----------------------------------------------------------%
% if visualize='true', user can define the time frame to visualize
t1=(datenum('20221201 000000','yyyymmdd HHMMSS')); % first timestep 
t2=(datenum('20221202 000000','yyyymmdd HHMMSS')); % last timestep 
dt=1/24;                      %delta t (day)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DeltaDir=360/nDir;
time=t1:dt:t2;
to=(datenum('1990-01-01 000000','yyyy-mm-dd HHMMSS',1990));
% gaugssian distribution for directional spread 
sigma=spread*DeltaDir; %degrees 
dir0=90:-DeltaDir:0;
dir0(end+1:nDir)=360-DeltaDir+dir0(end):-DeltaDir:90+DeltaDir-dir0(end);
dir=pi*dir0/180;
%----------------------------------------------------------%
% frequency (log)
Omega(1)=2*pi*f0;
for i=2:nfreq
Omega(i)=inc*Omega(i-1);
end
freq=Omega/2/pi;
%----- Jonswap spectrum -----------------------------------%
display('Generate Jonswap Spec ...')
g=9.81;
Gamma=3.3;   %peakedness parameter
Beta=5/4;
SigmaA=0.07;  %spectral width parameter
SigmaB=0.09;  %spectral width parameter
Omegap    = 2*pi/Tp;
sig = (Omega<=Omegap)*SigmaA+(Omega>Omegap)*SigmaB;
A     = exp(-((Omega/Omegap-1)./(sig*sqrt(2))).^2);  
alphabar = 5.058*(1-.287*log(Gamma))*(Hs/Tp^2)^2  ;                                     %modified Phillips constant
S = alphabar*g^2 .* Omega.^-5 .* exp(-(Beta*(Omega/Omegap).^-4)) .* Gamma.^A;      %spectra m^2.s
S = (2*pi)*S;
%------------------------------------------------------%
%check Hs
HScheck1=4.004*sqrt(abs(trapz(2*pi*freq,S)))/sqrt(2*pi);
%directional distribution
func=@(d)((1./(sigma*sqrt(2*pi))).*exp(-0.5*((d)/sigma).^2));
dd=-180+DeltaDir/2:DeltaDir:180-DeltaDir/2;
    for j=1:length(dd)
      A(j)=func(dd(j));
    end
    dd=dd+MeanDir;
    dd(dd>360)=dd(dd>360)-360;
    dd(dd<0)=dd(dd<0)+360;
 AA=interp1(dd,A,dir0,'linear','extrap');
    
for i=1:length(S)
      SPEC(:,i)=S(i).*AA;
end

for i=1:length(freq)
    Q(1,i) = abs(trapz(0:pi*DeltaDir/180:2*pi-pi*DeltaDir/180,SPEC(:,i)));
end
coef=nanmean(Q.\S);
Q=coef*Q;
HScheck2=4.004*sqrt(trapz(2*pi*freq,Q))/sqrt(2*pi);
SPEC=coef*SPEC;
%----------------------------------------------------------%
%repeat spectrum
for i=1:length(time)
    EFTH(:,:,1,i)=SPEC;
end
%----------------------------------------------------------%
% The boundary requres coordinates (lat, lon) depth, wind speed
% wind direction, current speed, current direction for each time step.
% Here wind and current are 0.
pointID = [pointID,repmat(' ', [1, 16-strlength(pointID)])];
Lon=Lon*ones(1,length(time));
Lat=Lat*ones(1,length(time));
curspd=curspd*ones(1,length(time));
curdir=curdir*ones(1,length(time));
wndspd=wndspd*ones(1,length(time));
wnddir=wnddir*ones(1,length(time));
dpt=dpt*ones(1,length(time));
%----------------------------------------------------------%
display (['Generating ', filename,' ...'])
%dump into netcdf
[filename] = write_directional_spectra_nc(filename,testcase,...
    pointID,Lat,Lon,dpt,wndspd,wnddir,curspd,curdir,time,freq,dir0,EFTH,coordinate);
%----------------------------------------------------------%         
tf=strcmp(visualize,'true');                 
if tf==1
display (['Plot Generation ', [filename,'.gif'],' ...'])
width=1280;  % Width of figure for movie [pixels]
height=880;  % Height of figure of movie [pixels]
left=200;     % Left margin between figure and screen edge [pixels]
bottom=200;  % Bottom margin between figure and screen edge [pixels]

h=figure('Color','w');
set(gcf,'Position', [left bottom width height])
for k=1:1:length(time)
    clf

sss1=subplot(2,1,1);

plot(freq,S,'-k','linewidth',2)
hold on
%plot(freq,Q,'--b','linewidth',2)
%hold on
xlabel('frequency','FontSize',12)
ylabel('mo','FontSize',12)

title(['H_s = ',num2str(Hs),' m; T_p = ', num2str(Tp),' s; Mean Dir = ',num2str(MeanDir),' deg; Directional Spread = ',num2str(sigma),' deg'])
box on
axis on
grid on
xlim([freq(1) freq(end)])

sss2=subplot(2,1,2);
SPECC=EFTH(:,:,1,k);
[~,cc] = polarPcolor(freq,[180*dir/pi 180*dir(1)/pi],...
     [SPECC; SPECC(1,:)],'Nspokes',36,'ncolor',10,'labelR','f (Hz)','Rscale','log');
 ylabel(cc,['Date = ',datestr(time(k))],'FontSize',14);

set(sss1, 'Position', [.08 .55 .88 .4]);
set(sss2, 'Position', [.08 .08 .88 .4]);
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
