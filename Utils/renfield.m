function stats = renfield(stats, oldname, newname)
%
%  stats = renfield(stats, fname)
%
% description: 
%    renames a field in the structure stats


if ~iscell(oldname)
   oldname = {oldname};
end
if ~iscell(newname)
   newname = {newname};
end

if length(oldname) ~= length(newname)
   error('renfield : incosnistent input!');
end

for i = 1:length(oldname)
   [stats.(newname{i})] = stats.(oldname{i});
   stats = rmfield(stats, oldname{i});
end

end