function [pol, varargout] = polygonFromMask(imgmsk, varargin)
%
% tri = polygonFromMask(img, varargin)
% 
% description: 
%      converts a image mask into a polygon / set of polygons
%
% input:
%      pol   cell array of countours not self-intersecting
%      param parameter struct with entries
%
% output:
%      pol   cell array of cell arrays representing the polygons
%
% See also: bwboundaries

%find boundaries
[b,~,~,tree] = bwboundaries(imgmsk);
pol = cellfunc(@(x) x', b)';

%tree if required
if nargout > 1
   varargout{1} = tree; 
end

end

