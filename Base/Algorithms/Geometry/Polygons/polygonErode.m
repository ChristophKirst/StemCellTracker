function [polbuf, varargout] = polygonErode(pol, dist, varargin) 
%
% pol = polygonBuffer(pol, varargin) 
%
% description: 
%     dilates or erodes the polygon
% 
% input:
%     pol     polygon as cell of oriented paths, each path is 2xn array of coords
%     dist    distance positive or negative
%     param   parameter struct with entries
%             .join     join type, 1=Square, 2=Round, 3=Miter (3)
%             .end      end type,  1=ClosedPolygon, 2=ClosedLine, 3=OpenButt, 4=OpenSquare, 5=OpenRound (1)
%             .scale    use this scale to convert the the plygon coords to integer ([]=automatic)
%             .simplify simplify polygon first to fix orientations for even odd interpretation (true)
%    
% output
%     polbuf  dilated polygon
%
% See also: polygonBuffer, polygonExecute

if nargout > 1
   [polbuf, varargout{1}] = polygonBuffer(pol, -dist, varargin{:});
else
   polbuf = polygonBuffer(pol, -dist, varargin{:});
end

end