function img = stitch2ImagesByOverwrite(img1, img2, shift, varargin)
%
% img = stitch2ImagesByOverwrite(img1, img2, shift)
%
% description:
%     stitches two images of same size together, completely overwriting img1 with img2 in overlap region
%
% input:
%     img1, img2   images to stitch together
%     shift        relative shift between img1 and img2 ([0,0] = no shift)
%
% output:
%     img          stitched image
%
% See also: stitchImages, alignImages

%transform relative shift into absolute shift for img1 and img2
dim = length(shift);
shift1 = zeros(1,dim);
shift2 = zeros(1,dim);

for d = 1:dim
   if shift(d) > 0
      shift2(d) = shift(d);
   else
      shift1(d) = -shift(d);
   end
end

siz = max([size(img1) + shift1; size(img2) + shift2]);
img = zeros(siz);
img = imreplace(img, img1, shift1 + 1);
img = imreplace(img, img2, shift2 + 1);

end

