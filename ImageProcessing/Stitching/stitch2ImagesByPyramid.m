function img = stitch2ImagesByPyramid(img1, img2, shift, varargin)
%
% img = stitch2ImagesByPyramid(img1, img2, shift)
%
% description:
%     stitches two images together using laplacian pyramid
%
% input:
%     img1, img2   images to stitch together
%     shift        relative shift between img1 and img2 ([0,0] = no shift)
%
% output:
%     img          stitched image
%
% See also: stitchImages, stitchImagesByMean, stitchImagesByWatershed, alignImages

img = stitchImagesByPyramid({img1, img2}, {zeros(1,ndims(shift)), shift}, varargin{:});

end

