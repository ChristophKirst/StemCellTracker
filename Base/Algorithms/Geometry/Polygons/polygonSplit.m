function [polys, varargout] = polygonSplit(poly, varargin) 
%
% polys = polygonSplit(pol, param) 
%
% description: 
%     splits the polygon pol into connected components
% 
% input:
%     pol     polygon as cell of oriented paths, each path is 2xn array of coords
%     param   parameter struct with entries
%             .filltype        EvenOdd = 1, NonZero = 2, Positive = 3, Negative = 4 (1)
%             .closed          default closed type of polygons (true)
%             .scale           use this scale to convert the the ploygon coords to integer ([]=automatic)
%             .full            split fully or only top level components (true)
%    
% output
%     polys   splitted polygons
%

[poly, tree] = polygonSimplify(poly, varargin{:});

polys = polygonTreeToCell(tree, varargin{:});

%tree if required
if nargout > 1
   varargout{1} = cellfunc(@(x) tree(x,x), polys); 
end

%transform numbers into polys
polys = cellfunc(@(x) poly(x), polys);

end
