function pol = polygonShift(pol, shift)
%
% pol = polygonShift(pol, shift)
% 
% description: 
%      shifts a polygon in space
%
% input:
%      pol   cell array of countours not self-intersecting
%      shift shift 
%
% output:
%      pol   cell array of cell arrays representing the polygons
%
% See also: polygonScale, polygonTransform

shift = shift(:);
pol = cellfunc(@(x) x + repmat(shift, 1, size(x,2)), pol);

end

