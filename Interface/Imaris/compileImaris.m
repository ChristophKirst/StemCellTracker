function compileImaris
%
% comileImaris()
%
% description:
%     compile files needed for the Imaris inteface
%

% create the statistics rgb color map
cmap = imarisrgbstatisticscolormap(8);
writeSSV(fullfile(fileparts(mfilename('fullpath')), 'colormap.pal'), floor(255 * cmap))

end