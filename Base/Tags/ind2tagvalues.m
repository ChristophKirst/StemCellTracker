function [tvs, ids] = ind2tagvalues(tagranges, i)

si = tagrangesize(tagranges);

if i > prod(si) || i < 1
   error('ind2tagvalues: index %g out of bounds %g:%g !', i, 1, prod(si));
end

ids = imind2sub(si, i);

names = fieldnames(tagranges);
nnames = length(names);

tvs = cell(1, nnames);

for i = 1:nnames
   tr = tagranges.(names{i});
   if iscell(tr)
      tvs{i} = tr{ids(i)};
   else
      tvs{i} = tr(ids(i));
   end
end

end