function img = stitchImagesByMax(imgs, shifts)
%
% img = stitchImagesByMax(imgs, shifts)
%
% description:
%     stitches images in cell array imgs together using shifts and max in overlapping 
%     regions
%
% input:
%     imgs         images as cell array
%     shifts       relative shifts between imgs{1} and imgs{k} ([0,0(,0)] = no shift)
%
% output:
%     img          stitched image
%
% See also: stitchImages, stitchImagesByMean, stitchImagesByMin, stitchImagesByOverwrite, alignImages

isizes = cellfunc(@size, imgs);
[ashifts, asize] = absoluteShiftsAndSize(shifts, isizes);

% find split into overlapping rectangles
[regs, ids] = findOverlapRegions(ashifts, isizes);


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

   img = imreplace(img, imgr, re(1,:));
   
end
      
end