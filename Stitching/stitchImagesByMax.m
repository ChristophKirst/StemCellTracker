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

imgsizes = cellfun(@size, imgs, 'UniformOutput', false);
[ashifts, asize] = absoluteShiftsAndSizes(shifts, imgsizes);


% calculate overlap regions -> simple but memory intensive
imgi = zeros(asize);
ind = 1;
for i = 1:numel(imgs)
   imgi = implus(imgi, ones(imgsizes{i}) * ind,  ashifts{i} + 1);
   ind = ind * 2;
end
inds = unique(imgi);
inds(inds == 0) = [];

% compose image
img = zeros(asize);

for i = inds'
   % overlapping region

   iovl = strfind(flip(dec2bin(i)), '1');
   posi = find(imgi == i);
   pos = imind2sub(asize, posi);
   

   % gather the overlap
   ovl = [];
   for j = 1:length(iovl)
      iimg = iovl(j);
      p0 = pos - repmat(ashifts{iimg}, size(pos, 1), 1);
      p0 = imsub2ind(imgsizes{iimg}, p0);
      ovl = cat(2, ovl, imgs{iimg}(p0));
   end 
   
   ovl = max(ovl, [], 2);
   img(posi) = ovl(:);
end

end