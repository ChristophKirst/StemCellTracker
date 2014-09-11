function tag = ind2tag(tagranges, i)


si = tagsize(tagranges);
ids = imind2sub(si, i);

names = fieldnames(tagranges);
nnames = length(names);

tvs = cell(1, 2* nnames);
for i = 1:nnames
   tvs{i} = names{i};
   tr = tagranges.(names{i});
   tvs{i+1} = tr{ids(i)};
end

tag = struct(tvs{:});

end