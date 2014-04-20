function name = tags2name(tfrmt, tags, tagnames)
%
% name = tags2name(tfrmt, tags, tagnames)
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
%    tags        struct/vector or cell with the image specifications
%    tagnames    (optional) list of names for the tag values and their ordering
%                if not specified this tagnames = tagformat2names(tfrmt)
%
% output:
%    name        the name with tag replaced by values
%
% See also: name2tags, tagformat, num2str0, tagformat2tagnames

name = tfrmt;

if nargin < 2 
   return 
elseif isempty(tags)
   if nargin == 2
      return
   elseif ~isempty(tagnames) 
      warning('tags2name: there are tagnames but no tag values'); 
      return
   end  
end

%convert tags to a struct that is a map string -> value
if ~isstruct(tags)
   if nargin == 2 || isempty(tagnames)
      tagnames = tagformat2tagnames(tfrmt);
   end

   if ~iscellstr(tagnames)
      error('tags2name: cannot infer tag names')
   end
   
   if ~iscell(tags)
      tags = num2cell(tags);
   end
   
   tags = {tagnames{:}; tags{:}};
   tags = struct(tags{:});
end


tagnames = fieldnames(tags);
for i = 1:length(tagnames)
   repl = regexp(name, ['<\s*?' tagnames{i} '\s*?>'],'match');
   if ~isempty(repl)
      name = strrep(name, repl{1}, num2str(tags.(tagnames{i})));
   end
   k = regexp(name, ['<\s*?' tagnames{i} '\s*?,\s*?(?<k>\d+)\s*?>'], 'names');
   if ~isempty(k)
      k = {k.k};
      for l = 1:length(k)
         repl = num2str0(tags.(tagnames{i}), str2double(k{l}));
         name = regexprep(name,  ['<\s*?' tagnames{i}  '\s*?,\s*?' k{l} '\s*?>'], repl);
      end
   end
end
     

end
   
   



