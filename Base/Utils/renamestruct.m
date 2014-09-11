function s = renamestruct(s, oldnames, newnames)
%
% s = renamestruct(s, oldnames, newnames)
%
% description:
%   renames fields with oldnames to newnames

if ischar(oldnames)
   oldnames = {oldnames};
end

if ischar(newnames)
   newnames = {newnames};
end

for i = 1:length(oldnames)
   %[s.(newnames{i})] = s.(oldnames{i});
   %s = rmfield(s,oldnames{i});
   
   s = renameStructField(s, oldnames{i}, newnames{i});
end

end