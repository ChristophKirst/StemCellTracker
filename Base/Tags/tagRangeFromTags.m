function tagRange = tagRangeFromTags(tags, varargin)
%
% tagranges = tagRangeFromTags(tags, varargin)
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
%      tagRange  struct with .name = {range} entries
%
% See also: tagExpression

param = parseParameter(varargin);

names = fieldnames(tags);
nnames = length(names);

tags(1)
tags(2)

tagRange = struct();
for i = 1:nnames
   tgs = {tags.(names{i})};
   
   if iscellstr(tgs)
      tagRange.(names{i}) = unique(tgs);
   else
      tagRange.(names{i}) = num2cell(unique([tags.(names{i})]));
   end
      
end


% simple but fast test for multiplicative data -> might miss certain very special cases
if getParameter(param, 'check', false)
   n = 1;
   for i = 1:nnames
      n = n * length(tagRange.(names{i}));
   end
   
   if length(tags) ~= n
      error('tagRangeFromTags: tags are not multiplicative!')
   end
end