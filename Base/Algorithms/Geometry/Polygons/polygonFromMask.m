function varargout = polygonFromMask(imgmsk, varargin)
%
% tri = polygonFromMask(img, varargin)
% 
% description: 
%      converts a image mask into a polygon / set of polygons
%
% input:
%      pol   cell array of countours not self-intersecting
%      param parameter struct with entries
%            .split    split into multiple polygons
%
% output:
%      pol   cell array of cell arrays representing the polygons
%
% See also: bwboundaries

varargout{:} = polygonFromLabeledImage(imgmsk, varargin{:});

varargout{1} = varargout{1}{1};

end

