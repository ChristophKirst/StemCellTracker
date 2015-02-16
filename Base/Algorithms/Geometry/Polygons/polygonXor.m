function [pol, varargout] = polygonXor(pol, polxor, varargin) 
%
% pol = polygonXor(pol, polclip) 
%
% description: 
%     calculates union between polygons
% 
% input:
%     pol     polygon as cell of oriented paths, each path is 2xn array of coords
%     polint  polygon to xor with 
%     param   parameter struct with entries as in polygonExecute
%
% output
%     pol     clipped polygon 
%

if nargout> 1
   [pol, varargout{1}] = polygonExecute(pol, polxor, varargin{:}, 'operator', 'Xor');
else
   pol = polygonExecute(pol, polxor, varargin{:}, 'operator', 'Xor');
end

end






