function [tag, ids] = ind2tag(tagranges, i)
%
% tag = ind2tag(tagranges, i)
%
% description:
%   converts tagranges and index i to the corresponding single tag

si = tagrangesize(tagranges);
ids = imind2sub(si, i);

names = fieldnames(tagranges);
nnames = length(names);

tvs = cell(1, 2* nnames);
for i = 1:nnames
   tvs{2*i-1} = names{i};
   tr = tagranges.(names{i});
   if iscell(tr)
      tvs{2*i} = tr{ids(i)};
   else
      tvs{2*i} = tr(ids(i));
   end
end

tag = struct(tvs{:});

end