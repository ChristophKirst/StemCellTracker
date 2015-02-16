function c = polygonTreeToCell(tree, varargin)
%
% c = polygonTreeToCell(tree, parents,varargin)
%
% description:
%     split tree matrix into a cell represenatation
%
% input:
%     tree    matrix represetning nesting structure of polys as returne de.g. by polygonSimplify
%     param   parameter struct with entries
%             .split     split nested polygons (true)
%
% output:
%     c       cell array of polygon ids that prepresent possible nested outlines


sf = getParameter(parseParameter(varargin), 'split', true);

n = length(tree);
[nid, pid]= find(tree);

parents = zeros(1,n);
parents(nid) = pid;

% top level enclosing components
ids = find(parents == 0);
c = num2cell(ids);

[c, ids] = addHoles(c, c, parents);

while any(~cellfun(@isempty, ids))
   [c, ids] = addPolys(c, ids, parents, sf);
   [c, ids] = addHoles(c, ids, parents);
end

%c = cellfunc(@(x) x', c);

end



function [polys, ids] = addHoles(polys, ids, parents)
   for i = 1:length(ids)
      ids{i} = find(ismember(parents, ids{i})); % find holes
      polys{i} = [polys{i}, ids{i}];            % add to polys
   end
end

function [polys, newids] = addPolys(polys, ids, parents, sf)
   newids = ids;
   for i = 1:length(ids)
      newpolys = find(ismember(parents, ids{i})); % find new polys
      
      if sf
         polys = [polys, num2cell(newpolys)]; %#ok<AGROW>
         newids = [newids, num2cell(newpolys)]; %#ok<AGROW>
         newids{i} = [];
      else
         polys{i} = [polys{i}, newpolys];
         newids{i} = newpolys;
      end

   end
end