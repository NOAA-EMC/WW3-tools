function [filename] = write_directional_spectra_ascii(filename,testcase,...
    pointID,Lat,Lon,dpt,wndspd,wnddir,curspd,curdir,time,freq,dir,EF)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function writes Directional Spectral density file in %
% WW3 ascii format                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Ali Abdolali April 2017 ali.abdolali@noaa.gov        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testcase name: i.e. Inlet test'
% time
% ntime: length of time
% nfreq: number of frequency band
% nDir: number of directional band
% pointnumber : number of points- use 1 for this variable
% freq: frequency (Hz)
% dir: direction (rad)
% pointID: the boundary point name: i.e. 'b42001' 
% Lat: latitude (degrees)
% Lon: longitude (degrees)
% dpt: depth (m) [pointnumber x 1]
% wndspd: wind  speed (m/s) [pointnumber, ntime]
% wnddir: wind direction (degrees) [pointnumber, ntime]
% curspd: current speed (m/s) [pointnumber, ntime]
% curdir: current direction (degrees) [pointnumber, ntime]
% EF: directional spectral density [nDir x nfreq  x pointnumber x ntime]
% filename: name of output file 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%example: [filename] = write_directional_spectra_ascii('B42001.spc',...
%'inlet','B42001',42,221.1,1521,30,270,0.5,36,time,freq,dir,EF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nDir,nfreq,pointnumber,ntime]=size(EF);


fileID = fopen([filename,'.spc'],'w');
fprintf(fileID,'%s     %d    %d     %d %s\n', 'WAVEWATCH III SPECTRA',...
                                         nfreq,nDir,pointnumber,testcase);

% write frequency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=floor((nfreq)/8);
remn=nfreq-(n*8);
for i=1:n
    fprintf(fileID,' %E',freq((i-1)*8+1:i*8));
        fprintf(fileID,'\n');
end
if remn>=1
  fprintf(fileID,' %E',freq(n*8+1:end));
          fprintf(fileID,'\n');
end

% write direction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   
n=floor((nDir)/7);
remn=nDir-(n*7);
for i=1:n
   fprintf(fileID,'  %E',dir((i-1)*7+1:i*7));
   fprintf(fileID,'\n');
end

if remn>=1
  fprintf(fileID,'  %E',dir(n*7+1:end));
          fprintf(fileID,'\n');
end

% write directional spectral density %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for it=1:ntime
    for ip=1:pointnumber
        

fprintf(fileID,'%s\n', datestr(time(it),formatOut));
fprintf(fileID,'%s  %2.2d     %2.2d   %3.3f %3.2d %1.2d %3.1d  %3.1d\n', ...
    pointID(ip),Lat(ip),Lon(ip),dpt(ip),wndspd(ip,it),wnddir(ip,it),...
                                          curspd(ip,it),curdir(ip,it));

 EFTH(:,:)=EF(:,:,ip,it);

      for i=1:nDir
          EFTHH(1,((i-1)*nfreq)+1:i*nfreq)=EFTH(i,:);
      end

      n=floor((nDir*nfreq)/7);
      remn=nDir*nfreq-(n*7);
      for i=1:n
        fprintf(fileID,'  %E',EFTHH((i-1)*7+1:i*7));
        fprintf(fileID,'\n');
      end

      if remn>=1
         fprintf(fileID,'  %E',EFTHH(n*7+1:end));
         fprintf(fileID,'\n');
      end
     end
end
fclose(fileID);
end