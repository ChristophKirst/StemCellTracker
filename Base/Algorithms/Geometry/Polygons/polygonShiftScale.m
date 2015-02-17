function pol = polygonShiftScale(pol, shift, scale)
%
% pol = polygonShiftScale(pol, shift)
% 
% description: 
%      shifts and then scales a polygon in space
%
% input:
%      pol   cell array of countours not self-intersecting
%      shift shift 
%      scale scale
%
% output:
%      pol   cell array of cell arrays representing the polygons
%
% See also: polygonScale, polygonTransform

shift = shift(:);
if length(scale) == 1
   scale = [scale, scale];
end
scale = scale(:);

pol = cellfunc(@(x) (x + repmat(shift, 1, size(x,2))) .* repmat(scale, 1, size(x,2)), pol);

end

