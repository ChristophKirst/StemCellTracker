function img = stitchImagesByMean(imgs, shifts, varargin)
%
% img = stitchImagesByMean(imgs, shifts)
%
% description:
%     stiches images in cell array imgs together using shifts and mean in overlapping  %     regions
%
% input:
%     imgs         images as cell array
%     shifts       relative shifts between imgs{1} and imgs{k} ([0,0(,0)] = no shift)
%
% output:
%     img          stiched image
%
% See also: stitchImages, stitchImagesByMax, stitchImagesByMin, stitchImagesByOverwrite, alignImages

isizes = cellfunc(@size, imgs);
[ashifts, asize] = absoluteShiftsAndSize(shifts, isizes);

% find split into overlapping rectangles
[regs, ids] = stitchImagesOverlapRegions(ashifts, isizes);

% compose image
img = zeros(asize);

for i = 1:length(regs)
   id = ids{i};
   re = regs{i};
   n = length(id);
 
   if n > 1
      si = re(2,:)-re(1,:)+1;
      imgr = zeros([n, si]);
      for k = 1:length(id)   
         sh = ashifts{id(k)};
         ie = imextract(imgs{id(k)}, re(1,:)-sh, re(2,:)-sh);
         imgr(k, :) = ie(:);
      end
      imgr = mean(imgr, 1);
      imgr = reshape(imgr(1,:), si);
   else
      sh = ashifts{id(1)};
      imgr = imextract(imgs{id(1)}, re(1,:)-sh, re(2,:)-sh);
   end
   
   img = imreplace(img, imgr, re(1,:));
   
end

end