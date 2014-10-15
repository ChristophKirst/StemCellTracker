function texpr = tagExpressionRename(texpr, oldnames, newnames)
%
% texpr = tagExpressionRename(texpr, oldnames, newnames)
%
% description:
%   renames tags with oldnames to newnames
%
% input:
%   texpr    tag expression
%   oldnames old names of tags in texpr as char or cellstr
%   newnames new names of tags in texpr as char or cellstr



if ischar(oldnames)
   oldnames = {oldnames};
end

if ischar(newnames)
   newnames = {newnames};
end


[tnames, tsplit, tinfo] = tagExpressionToTagNames(texpr);

[tids, nids] = ismember(tnames, oldnames);
tids = find(tids);

for i = 1:length(tids)
   ti = tids(i);
   tinfo(ti).tag = strrep(tinfo(ti).tag, tinfo(ti).name, newnames{nids(i)});
   tinfo(ti).name = newnames{nids(i)};
end

texpr = tagInfoToTagExpression(tsplit, tinfo);

end