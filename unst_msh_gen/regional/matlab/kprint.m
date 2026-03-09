function kprint(fl);
% function kprint(fl);
% prints and trims (using convert) current figure to file name fl
% i.e. kprint('MyGraphicFile.jpg')
h=gcf;
saveas(h, fl);
system(['convert -trim ',fl,' ',fl]);
