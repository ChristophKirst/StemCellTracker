function [polbuf, varargout] = polygonBuffer(pol, dist, varargin) 
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
%     polbuf  bufered polygon
%

param = parseParameter(varargin);

joinTypes = {'Square', 'Round', 'Miter'};
endTypes  = {'ClosedPolygon', 'ClosedLine', 'OpenButt', 'OpenSquare', 'OpenRound'};

[~,jT] = getParameterInNames(param, 'join', joinTypes, 3);
[~,eT] = getParameterInNames(param, 'end' , endTypes, 1);
sc = getParameter(param, 'scale', []);

adist = abs(dist);

if ~iscell(pol)
   pol = {pol};
end

for i = 1:length(pol)
   if size(pol{i},1) ~= 2
      error('polygonBuffer: polygon coordinate size %d is not 2 in contour %d!', size(pol{i},1), i);
   end
end

if getParameter(param, 'simplify', true)
   pol = polygonSimplify(pol);
end

if isempty(sc)
   % check if we have integer array
   if all(cellfun(@isinteger, pol)) % || all(cellfun(@(x) all(rem(x(:),1)==0), pol))
      sc = 1;
   else
      pmax = max(cellfun(@(x) max(x(:)), pol)) + adist;
      pmin = min(cellfun(@(x) max(x(:)), pol)) - adist;
      
      sc =  9223372036854775807 / 1000 / (pmax-pmin + 1);
   end
end

polbuf = cellfunc(@(x) int64(sc*x), pol);

if nargout > 1
   [polbuf, varargout{1}] = mexPolygonBuffer(polbuf, dist * sc, jT-1, eT-1);
else
   polbuf = mexPolygonBuffer(polbuf, dist * sc, jT-1, eT-1);
end

polbuf = cellfunc(@(x) double(x)/sc, polbuf);

end






