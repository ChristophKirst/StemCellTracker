function tagRange = tagRangeFromFirstAndLastString(tagExpr, strFirst, strLast)
%
% tagranges = tagRangeFromFirstAndLastString(tags, strFirst, strLast)
%
% description:
%      converts two strings to tag ranges assuming numeric values for all tags and increments of 1
%
% input:
%      tagExpr   the tag expression
%      str*      first and last strings
%      param     additional tag range specs
%
% output:
%      tagRange  struct with .name = {range} entries
%
% See also: tagExpression

tagFirst = tagExpressionToTags(tagExpr, strFirst);
tagLast  = tagExpressionToTags(tagExpr, strLast);

tnf = fieldnames(tagFirst);
tnl = fieldnames(tagLast);

if ~isequal(tnf, tnl)
   error('tagRangeFromFirstAndLastString: inconsistent first and last strings %s and %s', strFirst, strLast);
end

for i = 1:length(tnf);
   fn = tnf{i};
   % error when these are strings
   vf = tagFirst.(fn);
   vl = tagLast.(fn);
   
   if ischar(vf)
      if ~ischar(vl)
         error('tagRangeFromFirstAndLastString: inconsistent classes for tag %s', fn);
      end
      
      tagRange.(fn) = {vf, vl};
   elseif isnumeric(vf)
      if ~isnumeric(vl)
         error('tagRangeFromFirstAndLastString: inconsistent classes for tag %s', fn);
      end
      
      tagRange.(fn) = num2cell((tagFirst.(fn)):(tagLast.(fn)));
   else
      error('tagRangeFromFirstAndLastString: inconsistent classes for tag %s', fn);
   end
end

end