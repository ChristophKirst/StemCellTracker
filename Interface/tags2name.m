function name = tags2name(tfrmt, tags)
%
% name = tags2name(tfrmt, tags)
%
% description:
%    generates a name from a tagged format string tfrmt and tag values tags
%    replacements are done as follows:
%    for a field with name xxx in the struct tags the substring 
%       <xxx>    is replaced by the value, 
%       <xxx,k>  is replaced by the value xxx using k digits with trailling zeros

%
% input:
%    tfrmt       the tagged string
%    tags        struct with the image specifications as tags.name -> val
%
% output:
%    name        the name with tag replaced by values
%
% See also: name2tags, tagformat, num2str0, tagformat2tagnames

name = tfrmt;

if nargin < 2 || isempty(tags) && isemptystruct(tags)
   return 
end

[tnames, tw, ~, torig] = tagformat2tagnames(tfrmt);

tagnames = fieldnames(tags);
name = tfrmt;

for i = 1:length(tagnames)
   k = find(ismember(tnames, tagnames{i}));
   if isempty(k)
      warning('tags2name: tag %s is not in the tag format', tagnames{i})
   else
      if tw{k} == 0
         name = strrep(name, torig{k}, num2str(tags.(tagnames{i})));
      else
         name = strrep(name, torig{k}, num2str0(tags.(tagnames{i}), tw{i}));
      end
   end
end   

end
   
   



