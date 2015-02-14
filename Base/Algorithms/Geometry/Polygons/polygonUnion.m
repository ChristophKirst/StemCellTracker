function polbuf = polygonUnion(pol, varargin) 
%
% pol = polygonBuffer(pol) 
%
% description: 
%     unions the polygons in pol
% 
% input:
%     pol     polygon as cell of oriented paths, each path is 2xn array of coords
%     param   parameter struct with entries
%             .filltype EvenOdd = 0, NonZero = 1, Positive = 2, Negative = 3 (0)
%             .scale   the routines convert the plygon to integer, use this scale to convert ([]=automatic)
%    
% output
%     polbuf  bufered polygon
%

param = parseParameter(varargin);

fT = getParameter(param, 'filltype', 0);
sc = getParameter(param, 'scale', []);

if ~iscell(pol)
   pol = {pol};
end

for i = 1:length(pol)
   if size(pol{i},1) ~= 2
      error('polygonBuffer: polygon coordinate size %d is not 2 in contour %d!', size(pol{i},1), i);
   end
end

if isempty(sc)
   pmax = max(cellfun(@(x) max(x(:)), pol));
   pmin = min(cellfun(@(x) max(x(:)), pol));
   
   sc =  9223372036854775807 / 1000 / (pmax-pmin + 1);
end

polbuf = cellfunc(@(x) int64(sc*x), pol);

polbuf = mexPolygonUnion(polbuf, fT, fT);

polbuf = cellfunc(@(x) double(x)/sc, polbuf);

end






