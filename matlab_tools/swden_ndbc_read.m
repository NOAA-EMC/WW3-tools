function [SWDEN] = swden_ndbc_read(ncf,deltatheta,theta0,freq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function reads NDBC directional spectral density file%
% netcdf format from https://dods.ndbc.noaa.gov/            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Ali Abdolali Feb 2023 ali.abdolali@noaa.gov          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  input data %--------------------------------------------%
% ncf: name of netcdf file
% deltatheta: direction resolution (degree)
% theta0: first dir (degree)
% freq: Frequency for interpolation (Hz)
%  output data %--------------------------------------------%
%There are two outputs: 
%1) Orig: on original frequency
%2) Int: Interpolated on input frequency

% time (Matlab time)
% f: frequency (Hz)
% Dir: direction (degree)
% DENS: Directional Spectral Density [freq,direction,stations,time]
% SPEC: Spectral Density [freq,stations,time]
% Hs: Significant Wave Heigth (m)
% Fp: Peak Freq (Hz)
%----------------------------------------------------------%
time=convert_time(ncf,'time');
f(1,:)=double(ncread(ncf,'frequency'));
r1=ncread(ncf,'wave_spectrum_r1');
R1(:,:)=r1(1,1,:,:);
r2=ncread(ncf,'wave_spectrum_r2');
R2(:,:)=r2(1,1,:,:);
spec=ncread(ncf,'spectral_wave_density');
SPEC(:,:)=spec(1,1,:,:);
alpha1=ncread(ncf,'mean_wave_dir');
a1(:,:)=alpha1(1,1,:,:);
alpha2=ncread(ncf,'principal_wave_dir');
a2(:,:)=alpha2(1,1,:,:);
theta(1,:)=[90-theta0:-deltatheta:0+theta0 360+theta0-deltatheta:-deltatheta:90+deltatheta-theta0];
Hs(:,1)=4.004*sqrt(abs(trapz(2*pi*f,SPEC(:,:))))/sqrt(2*pi);

SPECmax=nanmax(SPEC);
for i=1:length(time)
    clear i1
    clear j1
    [i1,j1]=find(SPEC(:,i)==SPECmax(i));
    Fp(i,1)=f(i1(1));
end

    
for i=1:length(theta)
  DENS(i,:,:)=SPEC.*(1/pi).*(0.5+R1.*cosd(theta(i)-180-a1)+R2.*cosd(2*(theta(i)-180-a2)));
end
%----------------------------------------------------------%
SWDEN.Orig.f=f;
SWDEN.Orig.SPEC=SPEC;
SWDEN.Orig.Dir=theta;
SWDEN.Orig.Hs=Hs;
SWDEN.Orig.Fp=Fp;
SWDEN.Orig.DENS=DENS;
SWDEN.Orig.time=time;
%----------------------------------------------------------%
%interpolate
R1Int=interp1(f,R1,freq);
R2Int=interp1(f,R2,freq);
a1Int=interp1(f,a1,freq);
a2Int=interp1(f,a2,freq);
SPECInt=interp1(f,SPEC,freq);
SPECInt(isnan(SPECInt))=0;
HsInt(:,1)=4.004*sqrt(abs(trapz(2*pi*freq,SPECInt(:,:))))/sqrt(2*pi);
clear SPECmax
SPECmax=nanmax(SPECInt);
for i=1:length(time)
    clear i1
    clear j1
    [i1,j1]=find(SPECInt(:,i)==SPECmax(i));
    FpInt(i,1)=freq(i1(1));
end

   
for i=1:length(theta)
  DENSInt(i,:,:)=SPECInt.*(1/pi).*(0.5+R1Int.*cosd(theta(i)-180-a1Int)+R2Int.*cosd(2*(theta(i)-180-a2Int)));
end
%----------------------------------------------------------%
DENSInt(isnan(DENSInt))=0;
SWDEN.Int.f=freq;
SWDEN.Int.SPEC=SPECInt;
SWDEN.Int.Dir=theta;
SWDEN.Int.Hs=HsInt;
SWDEN.Int.Fp=FpInt;
SWDEN.Int.DENS=DENSInt;
SWDEN.Int.time=time;
%----------------------------------------------------------%