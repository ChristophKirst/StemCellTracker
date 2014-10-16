function img = stitchImagesByOverwrite(imgs, ipos, varargin)
%
% img = stitchImagesByOverwrite(imgs, shifts)
%
% description:
%     stitches images together completely overwriting earlier images in imgs cell array
%
% input:
%     imgs         images to stitch together in cell array
%     param        parameter struct with entries
%                  .size     final image size, ipos are interpreted as shifts if s other wise as abaolute positions in image of this size ([])
%
% output:
%     img          stitched image
%
% See also: stitchImagesByMean, stitchImagesByMax, stitchImagesByMin, alignImages

param = parseParameter(varargin);

asize = getParameter(param, 'size', []);

if isempty(asize)
   isizes = cellfun(@size, imgs, 'UniformOutput', false);
   [ashifts, asize] = absoluteShiftsAndSize(ipos, isizes);
else
   ashifts = ipos;
end

img = zeros(asize);
for i = 1:numel(imgs)
   img = imreplace(img, imgs{i},  ashifts{i} + 1, 'chop', true);
end

end

