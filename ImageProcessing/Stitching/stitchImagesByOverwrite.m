function img = stitchImagesByOverwrite(imgs, shifts)
%
% img = stitchImagesByOverwrite(imgs, shifts)
%
% description:
%     stitches images together completely overwriting earlier images in imgs cell array
%
% input:
%     imgs         images to stitch together in cell array
%     shifts       relative shifts between imgs{1} and imgs{k} ([0,0(,0)] = no shift)
%
% output:
%     img          stitched image
%
% See also: stitchImagesByMean, stitchImagesByMax, stitchImagesByMin, alignImages

imgsizes = cellfun(@size, imgs, 'UniformOutput', false);
[ashifts, asize] = absoluteShiftsAndSize(shifts, imgsizes);

img = zeros(asize);
for i = 1:numel(imgs)
   img = imreplace(img, imgs{i},  ashifts{i} + 1);
end

end

