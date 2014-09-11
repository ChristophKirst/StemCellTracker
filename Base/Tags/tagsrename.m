function tags = tagsrename(tags, oldnames, newnames)
%
% tags = tagsrename(tags, oldnames, newnames)
%
% description:
%   renames tags with oldnames to newnames

tags = renamestruct(tags, oldnames, newnames);

end