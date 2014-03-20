function relabel = imrelabel(label)
%
% relabel = imrelabel(label)
%
% description:
%    relabels the labeled image using labels 1 to number of labeled regoins
%
% input:
%    label    labeled image (2D/3D)
%
% output:
%    relabel  relabeld image
%
% See also: imlabelseparate

idx = regionprops(label, 'PixelIdxList');
idx = {idx.PixelIdxList};

relabel = label;
k = 1;
for i = 1:length(idx)
   if ~isempty(idx{i})
      relabel(idx{i}) = k;
      k = k + 1;
   end
end

end
   