function texpr = tagexprrename(texpr, oldnames, newnames)
%
% texpr = tagexprrename(texpr, oldnames, newnames)
%
% description:
%   renames tags with oldnames to newnames
%
% input:
%   texpr    tag expression
%   oldnames old names of tags in texpr as char or cellstr
%   newnames new names of tags in texpr as char or cellstr
%
% See also: tagexpr2tagnames


if ischar(oldnames)
   oldnames = {oldnames};
end

if ischar(newnames)
   newnames = {newnames};
end


[tnames, tsplit, tinfo] = tagexpr2tagnames(texpr);

[tids, nids] = ismember(tnames, oldnames);
tids = find(tids);

for i = 1:length(tids)
   ti = tids(i);
   tinfo(ti).tag = strrep(tinfo(ti).tag, tinfo(ti).name, newnames{nids(i)});
   tinfo(ti).name = newnames{nids(i)};
end

texpr = taginfo2tagexpr(tsplit, tinfo);

end