function cmap = imarisrgbstatisticscolormap(nrgb)
%
% cmap = imarisrgbstatisticscolormap(nrgb)
%
% description:
%     create color table for usage with imarisrgb2statisticscolor
%
% input:
%     nrgb    (optional) number of steps in each color r,g,b (256)
%
% ouput:
%     cmap    a color map
%
% See also: imarisrgb2statistics

if nargin < 1
   nrgb = 8;
end

[r,g,b] = ndgrid(0:(nrgb-1), 0:(nrgb-1), 0:(nrgb-1));

cmap = [r(:), g(:), b(:)] / (nrgb-1);

end





