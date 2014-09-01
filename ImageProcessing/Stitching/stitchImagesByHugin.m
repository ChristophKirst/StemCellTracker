function img = stitchImagesByHugin(imgs, shifts, varargin)
%
% img = stitchImagesByHugin(imgs, shifts, param)
%
% description:
%     stiches images in cell array imgs together using Hugin image sticher (enblend)
%
% input:
%     imgs         images as cell array
%     shifts       relative shifts between imgs{1} and imgs{k} ([0,0(,0)] = no shift)
%     param        parameter as in histitch
%
% output:
%     img          stitched image
%
% See also: stitchImages, alignImages, histitch

img = histitch(imgs, shifts, varargin{:});

end