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
   
   s = renameSfield(s, oldnames{i}, newnames{i});
end

end


%matlabs function renameStructField is buggz!
function str = renameSfield(str, oldFieldName, newFieldName)

   if ~strcmp(oldFieldName, newFieldName)
      allNames = fieldnames(str);
      % Is the user renaming one field to be the name of another field?
      % Remember this.
      isOverwriting = ~isempty(find(strcmp(allNames, newFieldName), 1));
      matchingIndex = find(strcmp(allNames, oldFieldName));
      if ~isempty(matchingIndex)
         allNames{matchingIndex(1)} = newFieldName;
         [str.(newFieldName)] = str.(oldFieldName);
         str = rmfield(str, oldFieldName);
         if (~isOverwriting)
            % Do not attempt to reorder if we've reduced the number
            % of fields.  Bad things will result.  Let it go.
            str = orderfields(str, allNames);
         end
      end
   end

end
