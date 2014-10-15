function fns = tagExpressionToFiles(texpr, varargin)
%
% fns = tagExpressionToFiles(texpr, param)
%
% description:
%      returns filnames that are in accordance with the tag format
%
% input:
%      texpr   tag expression
%      param   (optional) parameter struct with entries
%              .check    check for consistent filenames if multi tags (false)
%
% output:
%      fns     filnames that match texpr
%
% See also: tagExpression


param = parseParameter(varargin);

re = tagExpressionToRegularExpression(texpr);
fns = tagExpressionToFileExpression(texpr);

fns = dirr(fns);

% sort out matching regexp
fns = fns(cellfun(@(x) ~isempty(regexp(x, re, 'once')), fns, 'UniformOutput', true));

if getParameter(param, 'check', false);
   %take care of multiple occurences of a tag
   [~,~,tinfo] = tagExpressionToTagNames(texpr);
   
   if isempty(fieldnames(tinfo))
      return
   end
   
   if any(cellfun(@length, {tinfo.pos}) > 1)
      tags = tagExpressionToTags(texpr, fns, 'check', false);
      fns2 = tagExpressionToString(texpr, tags);
      fns = fns(cellfun(@isequal, fns, fns2));
   end
end

end