function [bnd, tree, stats] = imlabelboundary(imglab, stats)
%
%  [imgout, stats] = imlabelboundary(imglab, stats)
%
% description:
%      returns bnd and tree for each label
%      as if bwboundaries was applied to each label dilated with [1 1; 1 1] 
%      dilation is used to make lines of pixels detectable by polygon conversion routines
%      boundaries are shifted by 0.5 to make polygon precisely border the pixels of lab
%
% input:
%      imglab    labeled image
%      stats     (optional) pre calculated image statistics
% 
% output:
%      bnd       cell array of cell arrays with boundaries for eahc label
%      tree      cell array of boundaries hierarchies
%      stats     calculated image statistics
%
% See also: imlabelbwboundaries, bwboundaries

isize = size(imglab);
dim = length(isize);

if dim ~= 2
   error('imlabelboundary: expect 2d labeled image');
end

if ~exist('stats', 'var')
   stats = imstatistics(imglab, {'PixelIdxList', 'BoundingBox'});
else
   stats = imstatistics(imglab, stats, {'PixelIdxList', 'BoundingBox'});
end

n = length(stats);

bnd = cell(1,n);
tree = cell(1,n);

imglab = padarray(imglab, [1,1]);

for l = 1:n
   if mod(l, 500) == 0
      fprintf('imlabelboundary: %g / %g\n', l, n);
   end
   
   idxpix = stats(l).PixelIdxList;
   
   if ~isempty(idxpix)
      ll = imglab(idxpix(1));
      
      
      stats(l).BoundingBox
      bbox = stats(l).BoundingBox
      bbox(3:4) = bbox(3:4) + 2

      obj = imextract(imglab, bbox);

      figure(6); clf; implot(obj)
      
      obj = imdilate(obj == ll, strel([1,1;1,1]), 'full'); % dilate each label to get lines etc for polygons right

      figure(7); clf; implot(obj, 'color.limit', [0,1])
      max(obj(:))
      
      [b, ~,~, tree{l}] = bwboundaries(obj);
      b
      
      bnd{l} = cellfunc(@(x) x' + repmat(bbox(1:dim) - 1.5, 1, size(x,1)), b'); %-1.5 to match polygon outline with pixel
   end
end

end