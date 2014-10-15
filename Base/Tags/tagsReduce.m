function [tags, texpr] = tagsReduce(tags, varargin)
%
% tags = tagsreduce(tags)
% [tags, texpr] = tagsReduce(tags, texpr)
%
% description:
%   checks for for identical tag values / tag ranges and only keeps the first
%   if texpr is given also removes tags not in the tag expression

names = fieldnames(tags)';

if nargin > 1
   texpr = varargin{1};
   [tnames, tsplit, tinfo] = tagExpressionToTagNames(texpr);
   
   % remove all field names in tags not in texpr
   ids = ismember(names, tnames);
   tags = rmfield(tags, names(~ids));
   names = fieldnames(tags);
   nnames = length(names);

   conn = zeros(nnames);
   for n1 = 1:nnames
      for n2 = n1+1:nnames 
         if isequal({tags.(names{n1})}, {tags.(names{n2})})
            conn(n1, n2) = 1;
         end
      end
   end
   comp = graphAdjacencyMatrixToConnectedComponents(conn);
   
   for i = length(comp):-1:1  
      taggroup  = tinfo(comp{i});
      
      tagsinfonew(i).name = taggroup(1).name;
      for l = 1:length(taggroup)
         taggroup(l).tag = strrep(taggroup(l).tag, taggroup(l).name, taggroup(1).name);
      end
      tagsinfonew(i).tag  = [taggroup.tag];
      tagsinfonew(i).pos  = [taggroup.pos];
      tagsinfonew(i).widht = [taggroup.width];
      tagsinfonew(i).type  = taggroup(1).type;
   end

   texpr = tagInfoToTagExpression(tsplit, tagsinfonew);
   
   ids = ismember(names, {tagsinfonew.name});
   tags = rmfield(tags, names(~ids));

else
   nnames = length(names);
   red = {};
   for n1 = 1:nnames
      for n2 = n1+1:nnames
         if isequal({tags.(names{n1})}, {tags.(names{n2})})
            red = [red, names(n2)]; %#ok<AGROW>
         end
      end
   end

   red = unique(red);
   tags = rmfield(tags, red); 
end

end