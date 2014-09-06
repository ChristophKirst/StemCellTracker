function [relabel, newstats] = imrelabel(imglab, stats)
%
% relabel = imrelabel(label)
%
% description:
%    relabels the labeled image using labels 1 to number of labeled regoins
%
% input:
%    imglab   labeled image (2D/3D)
%
% output:
%    relabel  relabeld image
%
% See also: imlabelseparate

if nargin < 2
   stats = regionprops(imglab, 'PixelIdxList');
else
   stats = imstatistics(imglab, stats, 'PixelIdxList');
end

idx = {stats.PixelIdxList};

relabel = imglab;
k = 1;
for i = 1:length(idx)
   if ~isempty(idx{i})
      relabel(idx{i}) = k;
      k = k + 1;
   end
end

if nargout > 1
   newstats = stats;
end

end
   