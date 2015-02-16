function [pol, varargout] = polygonIntersection(pol, polint, varargin) 
%
% pol = polygonIntersection(pol, polclip) 
%
% description: 
%     calculates union between polygons
% 
% input:
%     pol     polygon as cell of oriented paths, each path is 2xn array of coords
%     polint  polygon to intersect with 
%     param   parameter struct with entries as in polygonExecute
%
% output
%     pol     clipped polygon 
%

if nargout> 1
   [pol, varargout{1}] = polygonExecute(pol, polint, varargin{:}, 'operator', 'Intersection');
else
   pol = polygonExecute(pol, polint, varargin{:}, 'operator', 'Intersection');
end

end






