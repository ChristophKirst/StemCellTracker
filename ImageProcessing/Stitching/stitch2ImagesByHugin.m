function img = stitch2ImagesByHugin(img1, img2, shift, varargin)
%
% img = stitch2ImagesByHugin(img1, img2, shift, param)
%
% description:
%     stiches two images together using Hugin image sticher (enblend)
%
% input:
%     img1, img2   images
%     shift        relative shift between img1 and img1 ([0,0(,0)] = no shift)
%     param        parameter as in histitch
%
% output:
%     img          stitched image
%
% See also: stitchImages, alignImages, histitch

img = histitch({img1, img2}, {zeros(1,ndims(shift)), shift}, varargin{:});

end