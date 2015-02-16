function [bnd, tree, stats] = imlabelbwboundary(imglab, stats)
%
%  [imgout, stats] = imlabelbwboundary(imglab, stats)
%
% description:
%      returns bnd and tree for each label
%      as if bwboundaries was applied on each label separately
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
% See also: bwboundaries

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

for l = 1:n
   if mod(l, 500) == 0
      fprintf('imlabelboundary: %g / %g\n', l, n);
   end
   
   idxpix = stats(l).PixelIdxList;
   
   if ~isempty(idxpix)
      ll = imglab(idxpix(1));
      
      bbox = stats(l).BoundingBox;
      
      %flat objects with trailing size 1 do not work with matlab !
      if bbox(dim)==bbox(2*dim)
         if bbox(2*dim) < isize(dim)
            bbox(2*dim) = bbox(dim) + 1;
         else
            bbox(dim) = bbox(dim) - 1;
            if bbox(dim) < 1
               bbox(dim) = 1;
            end
         end
      end

      obj = imextract(imglab, bbox);
      %obj = imdilate(obj == ll, strel([1,1;1,1])); % dilate each label to get lines etc for polygons right

      [b, ~,~, tree{l}] = bwboundaries(obj);
      bnd{l} = cellfunc(@(x) x' + repmat(bbox(1:dim)-1,1, size(x,1)), b);
      
   end
end

end