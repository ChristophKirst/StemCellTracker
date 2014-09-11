function tagranges = tags2tagranges(tags, varargin)
%
% tagranges = tags2tagranges(tags, varargin)
%
% description:
%      converts tags to tag ranges
%
% input:
%      tags    struct array of tags
%      param   (optional) parameter struct with entries
%              .check    check for multiplicative tag grid (false)
%
% output:
%      tagranges  struct with .name = {range} entries
%
% See also: tagformat

param = parseParameter(varargin);

names = fieldnames(tags);
nnames = length(names);

tagranges = struct();
for i = 1:nnames
   tgs = {tags.(names{i})};
   
   if iscellstr(tgs)
      tagranges.(names{i}) = unique(tgs);
   else
      tagranges.(names{i}) = num2cell(unique([tags.(names{i})]));
   end
      
end


% simple but fast test for multiplicative data -> might miss certain very special cases
if getParameter(param, 'check', false)
   n = 1;
   for i = 1:nnames
      n = n * length(tagranges.(names{i}));
   end
   
   if length(tags) ~= n
      error('tags2tagranges: tags are not multiplicative!')
   end
end