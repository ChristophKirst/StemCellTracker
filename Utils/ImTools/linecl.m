function lc = linecl(x, y, varargin)
%
% lc = linecl(x,y,varargin)
%
% description:
%     creates a closed line
%
% input:
%    x,y      x,z coordinates 
%    varargin args passed to line
%
% output:
%    lc       closed line
%
% See also: line

   lc = line([x, x(1)], [y y(1)], varargin{:});
   
end