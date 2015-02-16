function [pol, varargout] = polygonDifference(pol, poldif, varargin) 
%
% pol = polygonDifference(pol, polclip) 
%
% description: 
%     calculates difference between polygons
% 
% input:
%     pol     polygon as cell of oriented paths, each path is 2xn array of coords
%     polclip the polygon used to define the clipping mask
%     param   parameter struct with entries as in polygonExecute
%
% output
%     polres  clipped polygon 
%

if nargout> 1
   [pol, varargout{1}] = polygonExecute(pol, poldif, varargin{:}, 'operator', 'Difference');
else
   pol = polygonExecute(pol, poldif, varargin{:}, 'operator', 'Difference');
end

end






