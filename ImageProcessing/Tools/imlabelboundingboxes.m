function [bboxes, varargout] = imlabelboundingboxes(label)
%
% bboxes = imlabelboundingboxes(label)
% [bboxeslo, bboxeshi] = imlabelboundingboxes(label)
%
% description:
%     finds the bouding box of each labeled segment in pixel coordinates
%
% input:
%     label       labeled image
% 
% output:
%     bboxes      bounding box as array of [bboxlo bboxhi]
%     bboxeslo    bounding box lower corners as array of [p,q,l]
%     bboxeshi    bounding box upper corner as array of [p,q,l]
%
% See also: imlabelseparate, regionprops

%%% regionprops version (faster)

dim = ndims(label);

bboxes = regionprops(label, 'BoundingBox');
bboxes = reshape([bboxes.BoundingBox], dim * 2, []);

ee = [2 1 3]; ee = ee(1:dim);
bboxes(1:dim, :) = round(bboxes(ee,:) + 0.5); % correct qp to pl and shift to full pixel

ee = ee + dim;
bboxeshi = bboxes(1:dim, :) + bboxes(ee, :) - 1;

if nargout  <= 1
   bboxes = vertcat(bboxes(1:dim, :), bboxeshi)';
else
   bboxes = bboxes(1:dim, :)';
   varargout{1} = bboxeshi';
end


%%% Version without regionprops - slow
% lab = imlabel(label);
% nlab = length(lab);
% dim = ndims(label);
% 
% bboxes = zeros(nlab,2* dim);
% i = 1;
% for l = lab
%    [mi, mx] = imboundingbox(label == l);
%    bboxes(i, :) = [mi, mx];
%    i = i + 1;
% end
% 
% if nargout > 1
%    varargout{1} = bboxes(:, dim+1:2*dim);
%    bboxes =  bboxes(:, 1:dim);
% end

end
