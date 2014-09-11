function tvs = ind2tagvalues(tagranges, i)

si = tagrangesize(tagranges);
ids = imind2sub(si, i);

names = fieldnames(tagranges);
nnames = length(names);

tvs = cell(1, nnames);

for i = 1:nnames
   tr = tagranges.(names{i});
   tvs{i} = tr{ids(i)};
end

end