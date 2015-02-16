function [pol, varargout] = polygonFromLabeledImage(imglab, varargin)
%
% tri = polygonFromLabeledImage(img, varargin)
% 
% description: 
%      converts a labeled image into a polygon / set of polygons
%
% input:
%      pol   cell array of countours not self-intersecting
%      param parameter struct with entries
%            .split   keep nested polygons or split (true)
%
% output:
%      pol   cell array of cell arrays representing the polygons
%
% See also: bwboundaries

%find boundaries
[b,tree] = imlabelboundary(imglab);

outree = {};
pol = {};
for i = 1:length(b)

   % hole structure
   polys = polygonTreeToCell(tree{i}, varargin{:});
   
   %tree if required
   if nargout > 1
      outree = [outree, cellfunc(@(x) tree(x,x), polys)];  %#ok<AGROW>
   end

   %transform numbers into polys
   pol = [pol, cellfunc(@(x) b{i}(x), polys)]; %#ok<AGROW>

end

%tree if required
if nargout > 1
   varargout{1} = outree;
end

end
