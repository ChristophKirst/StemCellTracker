function tagRange = tagRangeFromTagExpression(tagExpr, varargin)
%
% tagranges = tagRangeFromFileExpression(tagExpr, tag range specs)
% tagranges = tagRangeFromFileExpression(tagExpr, names, tag range specs)
%
% description:
%      converts tagExpr to tag range using thenames or if not given the file names matching tagExpr
%
% input:
%      tags    struct array of tags
%      param   (optional) parameter struct with entries
%              .fast     use first and last name only to infer tag range
%              other ranges specifed by user
%
% output:
%      tagRange  struct with .name = {range} entries
%
% See also: tagExpression

if nargin > 1 && iscellstr(varargin)
   names = varargin{1};
   varargin = varargin(2:end);
else
   names = [];
end

tagRange = parseParameter(varargin);
fast = getParameter(tagRange, 'fast', true);
check = getParameter(tagRange, 'check', false);

rmf = {'fast', 'check'};
for i = 1:length(rmf)
   if isfield(tagRange, rmf{i})
      tagRange = rmfield(tagRange, rmf{i});
   end
end

if isempty(names)
   fileExpr = tagExpressionToFileExpression(tagExpr, tagRange);
   names = dirr(fileExpr);
   %length(names)
   
   if isempty(names)
      error('tagRangeFromTagExpression: no files of form %s !', fileExpr);
   end
   
end
 
if isempty(names)
    error('tagRangeFromTagExpression: no names found !');
end


% for names with non trailing zeros the first and last file might no reflect the range
ll = unique(cellfun(@length, names));

if fast && length(ll) == 1
   tagRange = parseParameter(tagRangeFromFirstAndLastString(tagExpr, names{1}, names{end}), tagRange);
else
   tags = tagExpressionToTags(tagExpr, names);
   tagRange = parseParameter(tagRangeFromTags(tags), tagRange);

   % simple but fast test for multiplicative data -> might miss certain very special cases
   if check
      n = 1;
      nnames = fieldnames(tagRange);
      for i = 1:nnames
         n = n * length(tagRange.(names{i}));
      end
      
      if length(tags) ~= n
         error('tagRangeFromTagExpression: tags are not multiplicative!')
      end
   end
   
end


end