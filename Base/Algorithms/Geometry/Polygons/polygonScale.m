function pol = polygonScale(pol, sc)
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

if length(sc) == 1
   sc = [sc,sc];
end
sc = sc(:);

pol = cellfunc(@(x) x .* repmat(sc, 1, size(x,2)), pol);

end

