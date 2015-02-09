function img = stitchImagesByInterpolate(imgs, ipos, varargin)
%
% img = stitchImagesByInterpolate(imgs, shifts)
%
% description:
%     stiches images in cell array imgs together using linear interpolation in overlaps  %     regions
%
% input:
%     imgs         images as cell array
%     ipos         positions or relative shifts between imgs{1} and imgs{k} ([0,0(,0)] = no shift)
%     param        parameter struct with entries
%                  .size     final image size, ipos are interpreted as shifts if no size is specified, otherwise as absolute positions in image of this size ([])
%
% output:
%     img          stiched image
%
% See also: stitchImages, stitchImagesByMax, stitchImagesByMin, stitchImagesByOverwrite, alignImages

if numel(imgs) == 1
   img = imgs{1};
   return
end

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

% find weights (distances form image centroids)
% speed up by using same weight if sizes are all the same
eq = cellfun(@(x) isequal(x, isizes{1}), isizes);
eq = all(eq(:));

if eq
   iids = 1;
else
   iids = 1:length(imgs);
end

wghts = cell(size(imgs));
for i = iids;
 
%    % find the center of the images
%    is = isizes{i};
%    nd = length(is);
%    cids = cell(1, nd);
%    for n = 1:nd
%       isn = is(n);
%       if mod(isn, 2) == 0
%          cids{n} = isn/2 + [0,1];
%       else
%          cids{n} = (isn+1)/2;
%       end
%    end
%    w = zeros(is);
%    w(cids{:}) = 1;
%    w = bwdist(w);
%    wghts{i} = w;

   % distances form border
   is = isizes{i};
   w = zeros(is);
   nd = length(is);
   for n = 1:nd;
      bids = repmat({':'}, 1, nd);
      bids{n} = [1, is(n)+1];
      w(bids{:}) = 1;
   end
   w = bwdist(w) + 1;
   %figure(21); clf; implot(w);
   wghts{i} = w; 
end

if eq
   wghts = repmat({w}, size(imgs));
end
      
% compose image
img = zeros(asize);

for i = 1:length(regs)
   id = ids{i};
   re = regs{i};
   n = length(id);
 
   if n > 1
      si = re(2,:)-re(1,:)+1;
      imgr = zeros([n, si]);
      imgw = zeros([n, si]);
      for k = 1:length(id)   
         sh = ashifts{id(k)};
         ie = imextract(imgs{id(k)}, re(1,:)-sh, re(2,:)-sh);
         iw = imextract(wghts{id(k)}, re(1,:)-sh, re(2,:)-sh);
         imgr(k, :) = ie(:);
         imgw(k, :) = iw(:);
      end
  
      imgws = sum(imgw, 1); 
      for k = 1:length(id)
         imgw(k,:) = imgw(k,:) ./ imgws(:)';
      end
      
      imgr = sum(imgr .* imgw, 1);
      imgr = reshape(imgr(1,:), si);
   else
      sh = ashifts{id(1)};
      imgr = imextract(imgs{id(1)}, re(1,:)-sh, re(2,:)-sh);
   end
   
   img = imreplace(img, imgr, re(1,:), 'chop', true);
   
end

end