function img = stitchImagesByMax(imgs, ipos, varargin)
%
% img = stitchImagesByMax(imgs, shifts)
%
% description:
%     stitches images in cell array imgs together using shifts and max in overlapping 
%     regions
%
% input:
%     imgs         images as cell array
%     ipos         positions or relative shifts between imgs{1} and imgs{k} ([0,0(,0)] = no shift)
%     param        parameter struct with entries
%                  .size     final image size, ipos are interpreted as shifts if s other wise as abaolute positions in image of this size ([])
%
% output:
%     img          stitched image
%
% See also: stitchImages, stitchImagesByMean, stitchImagesByMin, stitchImagesByOverwrite, alignImages

isizes = cellfunc(@size, imgs);

param = parseParameter(varargin);

asize = getParameter(param, 'size', []);

if isempty(asize)
   [ashifts, asize] = absoluteShiftsAndSize(ipos, isizes);
else
   ashifts = ipos;
end

% find split into overlapping rectangles
[regs, ids] = stitchImagesOverlapRegions(ashifts, isizes);

% compose image
img = zeros(asize);

for i = 1:length(regs)
   id = ids{i};
   re = regs{i};

   sh = ashifts{id(1)};
   
   imgr = imextract(imgs{id(1)}, re(1,:)-sh, re(2,:)-sh);  
   for k = 2:length(id)   
      sh = ashifts{id(k)};
      imgr = max(imgr, imextract(imgs{id(k)}, re(1,:)-sh, re(2,:)-sh));
   end

   img = imreplace(img, imgr, re(1,:), 'chop', true);
   
end
      
end