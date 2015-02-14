function polbuf = polygonBuffer(pol, dist, varargin) 
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
%             .join    join type, 0=Square, 1=Round, 2=Miter (2)
%             .end     end type,  0=ClosedPolygon, 1=ClosedLine, 2=OpenButt, 3=OpenSquare, 4=OpenRoun
%             .scale   the routines convert the plygon to integer, use this scale to convert ([]=automatic)
%    
% output
%     polbuf  bufered polygon
%

param = parseParameter(varargin);

jT = getParameter(param, 'join', 2);
eT = getParameter(param, 'end' , 0);
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

if isempty(sc)
   pmax = max(cellfun(@(x) max(x(:)), pol)) + adist;
   pmin = min(cellfun(@(x) max(x(:)), pol)) - adist;
   
   sc =  9223372036854775807 / 1000 / (pmax-pmin + 1);
end

polbuf = cellfunc(@(x) int64(sc*x), pol);

polbuf = mexPolygonBuffer(polbuf, dist * sc, jT, eT);

polbuf = cellfunc(@(x) double(x)/sc, polbuf);

end






