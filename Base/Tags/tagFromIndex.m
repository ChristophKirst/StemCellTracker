function [tag, ids] = tagFromIndex(tagRange, id)
%
% tag = tagFromIndex(tagRange, id)
%
% description:
%   converts index id to the corresponding tag in the tagRange

si = tagRangeSize(tagRange);

names = fieldnames(tagRange);
nnames = length(names);

ids = imind2sub(si, id);

for j = length(id):-1:1
   tvs = cell(1, 2* nnames);
   for i = 1:nnames
      tvs{2*i-1} = names{i};
      tr = tagRange.(names{i});
      if iscell(tr)
         tvs{2*i} = tr{ids(j,i)};
      else
         tvs{2*i} = tr(ids(j,i));
      end
   end

   tag(j) = struct(tvs{:});
end

end