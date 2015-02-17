function [pol, varargout] = polygonFromLabeledImage(imglab, varargin)
%
% [pol, ids, tree] = polygonFromLabeledImage(imglab, varargin)
% 
% description: 
%      converts a labeled image into a polygon / set of polygons
%
% input:
%      pol   cell array of countours not self-intersecting
%      param parameter struct with entries
%            .split   split different polygons for label with disconnected components (false)
%
% output:
%      pol   cell array of cell arrays representing the polygons
%      ids   labels the polygons are associated to
%      tree  tree structure of polygons
%
% See also: bwboundaries

%find boundaries

param = parseParameter(varargin);

sp = getParameter(param, 'split', false');

[b,tree] = imlabelbwboundary(imglab);

outtree = {};
pol = {};
ids = [];

for i = 1:length(b)

   if ~isempty(b{i})
      

      % splitting
      if sp
         polys = polygonTreeToCell(tree{i}, varargin{:});
         
         %tree if required
         if nargout > 2
            treei = tree{i};
            outtree = [outtree, cellfunc(@(x) treei(x,x), polys)];  %#ok<AGROW>
         end

      else
         polys = {1:length(b{i})};
         
         if nargout > 2
            outtree = [outtree, tree(i)]; %#ok<AGROW>
         end
      end
      
      % transform pixel boundaries into polygons
      b{i} = cellfunc(@(x) polygonBuffer(x, 0.5, 'simplify', false, 'end', 2), b{i}); % potential clipper bug on left hand pixels -> use end=2
      b{i} = cellfunc(@(x) x{1}, b{i}); % keep outer boundary only
      
      newpol = cellfunc(@(x) b{i}(x)', polys);
      pol = [pol, newpol]; %#ok<AGROW>
      
      if nargout > 1
         ids = [ids, i * ones(1,length(newpol))]; %#ok<AGROW>
      end
      
   end
end

if nargout > 1
   varargout{1} = ids;
end

if nargout > 2
   varargout{2} = outtree;
end

end
