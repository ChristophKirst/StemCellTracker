function [poly, varargout] = polygonExecute(poly, polyclip, varargin) 
%
% [pol, tree] = polygonExecute(pol, polclip, param) 
%
% description: 
%     clips the polygon pol to to regoin defined by the polygon polclip using the defined operation
% 
% input:
%     pol     polygon as cell of oriented paths, each path is 2xn array of coords
%     polclip the polygon used to define the clipping mask
%     param   parameter struct with entries
%             .filltype.subject   EvenOdd = 1, NonZero = 2, Positive = 3, Negative = 4 (1)
%             .filltype.clip      EvenOdd = 1, NonZero = 2, Positive = 3, Negative = 4 (1)
%             .operator           Intersection = 1, Union = 2, Difference = 3, Xor = 5 (2)
%             .closed             default closed type of polygons (true)
%             .scale              use this scale to convert the the ploygon coords to integer ([]=automatic)
%    
% output
%     pol     clipped polygon
%     tree    hierarchical structure of polygons as matrix A as e.g. in bwboundaries
%

param = parseParameter(varargin);

fillTypes = {'EvenOdd', 'NonZero', 'Positive', 'Negative'};

[~,fTsubj] = getParameterInNames(param, 'filltype.subject', fillTypes, 1);
[~,fTclip] = getParameterInNames(param, 'filltype.clip', fillTypes, 1);

opTypes = {'Intersection', 'Union', 'Difference', 'Xor'};
[~, op] = getParameterInNames(param, 'operator', opTypes, 2);

sc = getParameter(param, 'scale', []);

cl = getParameter(param, 'closed', true);


%handle Polygon class here to get closed property
if ~iscell(poly)
   poly = {poly};
end

for i = 1:length(poly)
   if size(poly{i},1) ~= 2
      error('polygonExecute: polygon coordinate size %d is not 2 in contour %d!', size(poly{i},1), i);
   end
end

if ~iscell(polyclip)
   polyclip = {polyclip};
end

for i = 1:length(polyclip)
   if size(polyclip{i},1) ~= 2
      error('polygonExecute: clipping polygon coordinate size %d is not 2 in contour %d!', size(polyclip{i},1), i);
   end
end

if isempty(sc)
   % check if we have integer array
   if all(cellfun(@isinteger, poly)) && ... % (all(cellfun(@isinteger, poly)) || all(cellfun(@(x) all(rem(x(:),1)==0), poly))) && ...
      all(cellfun(@isinteger, polyclip)) % (all(cellfun(@isinteger, polyclip)) || all(cellfun(@(x) all(rem(x(:),1)==0), polyclip)))  
      sc = 1;
   else
      pmax = max(cellfun(@(x) max(x(:)), poly));
      pmin = min(cellfun(@(x) max(x(:)), poly));
      
      sc =  9223372036854775807 / 1000 / (pmax-pmin + 1);
   end
end

poly = cellfunc(@(x) int64(sc*x), poly);
polyclip = cellfunc(@(x) int64(sc*x), polyclip);

if nargout == 1
   poly = mexPolygonExecute(op-1, [fTsubj, fTclip]-1, poly, cl, polyclip, cl);
else
   [poly, varargout{1}] = mexPolygonExecute(op-1, [fTsubj, fTclip]-1, poly, cl, polyclip, cl);
end

poly = cellfunc(@(x) double(x)/sc, poly);

end






