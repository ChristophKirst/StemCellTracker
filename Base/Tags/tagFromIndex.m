function [tag, ids] = tagFromIndex(tagRange, i)
%
% tag = tagFromIndex(tagRange, i)
%
% description:
%   converts index i to the corresponding tag in the tagRange

si = tagRangeSize(tagRange);
ids = imind2sub(si, i);

names = fieldnames(tagRange);
nnames = length(names);

tvs = cell(1, 2* nnames);
for i = 1:nnames
   tvs{2*i-1} = names{i};
   tr = tagRange.(names{i});
   if iscell(tr)
      tvs{2*i} = tr{ids(i)};
   else
      tvs{2*i} = tr(ids(i));
   end
end

tag = struct(tvs{:});

end