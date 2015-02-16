function [poly, varargout] = polygonSimplify(poly, varargin) 
%
% [pol, tree] = polygonSimplify(pol, param) 
%
% description: 
%     clips the polygon pol to to regoin defined  by the polygon polclip
% 
% input:
%     pol     polygon as cell of oriented paths, each path is 2xn array of coords
%     param   parameter struct with entries
%             .filltype        EvenOdd = 1, NonZero = 2, Positive = 3, Negative = 4 (1)
%             .closed          default closed type of polygons (true)
%             .scale           use this scale to convert the the ploygon coords to integer ([]=automatic)
%    
% output
%     pol     clipped polygon
%     tree    hierarchical structure of polygons as matrix A as e.g. in bwboundaries
%

param = parseParameter(varargin);

fillTypes = {'EvenOdd', 'NonZero', 'Positive', 'Negative'};
[~, fTsubj] = getParameterInNames(param, 'filltype.subject', fillTypes, 1);

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

if isempty(sc)
   % check if we have integer array
   if (all(cellfun(@isinteger, poly))) && ... % || all(cellfun(@(x) all(rem(x(:),1)==0), poly))) && ...
      (all(cellfun(@isinteger, polyclip)))    % || all(cellfun(@(x) all(rem(x(:),1)==0), polyclip)))  
      sc = 1;
   else
      pmax = max(cellfun(@(x) max(x(:)), poly));
      pmin = min(cellfun(@(x) max(x(:)), poly));
      
      sc =  9223372036854775807 / 1000 / (pmax-pmin + 1);
   end
end

poly = cellfunc(@(x) int64(sc*x), poly);

if nargout == 1
   poly = mexPolygonExecute(1, [fTsubj, fTsubj]-1, poly, cl);
else
   [poly, varargout{1}] = mexPolygonExecute(1, [fTsubj, fTsubj]-1, poly, cl);
end

poly = cellfunc(@(x) double(x)/sc, poly);

end






