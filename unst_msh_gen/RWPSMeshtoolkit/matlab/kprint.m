function kprint(fl);
%eval(['print -djpeg ',fl]);
h=gcf;
saveas(h, fl);
system(['convert -trim ',fl,' ',fl]);
