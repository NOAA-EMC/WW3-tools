function kprint(fl);
% function kprint(fl);
% prints and trims (using convert) current figure
h=gcf;
saveas(h, fl);
system(['convert -trim ',fl,' ',fl]);
