function img = stitch2ImagesByMax(img1, img2, shift, varargin)
%
% img = stitch2ImagesByMax(img1, img2, shift)
%
% description:
%     stitches two images of same size together, calculating 
%     the max in the overlapping area
%
% input:
%     img1, img2   images to stitch together
%     shift        relative shift between img1 and img2 ([0,0] = no shift)
%
% output:
%     img          stitched image
%
% See also: stitchImages, stitchImagesByMean, stitchImagesByOverwrite, alignImages

%transform relative shift into absolute shift for img1 and img2
dim = length(shift);


% find overlap region
si1 = size(img1);
si2 = size(img2);

ovl1 = cell(1,dim);
ovl2 = cell(1,dim);

shift1 = zeros(1,dim);
shift2 = zeros(1,dim);

for d = 1:dim
   sh = shift(d);
   if sh > 0
      ovl1{d} = (sh + 1) : min(si1(d), si2(d)+sh);
      ovl2{d} = 1:(ovl1{d}(end)-ovl1{d}(1) + 1);
      shift2(d) = shift(d);
   else
      ovl2{d} = (-sh + 1) : min(si2(d), si1(d)-sh);
      ovl1{d} = 1:(ovl2{d}(end)-ovl2{d}(1) + 1);
      shift1(d) = -shift(d);
   end 
end

iovl1 = img1(ovl1{:});
iovl2 = img2(ovl2{:});
img2(ovl2{:}) = max(iovl1,iovl2);

siz = max([size(img1) + shift1; size(img2) + shift2]);
img = zeros(siz);
img = imreplace(img, img1, shift1 + 1, 'chop', true);
img = imreplace(img, img2, shift2 + 1, 'chop', true);

end

