function fn = tagExpressionToFileExpression(texpr, varargin)
%
% fn = tagExpressionToFileExpression(texpr)
% fn = tagExpressionToFileExpression(texpr, tagspecs)
%
% description:
%    returns filename pattern that is in accordance with the tag format
%
% input:
%    texpr   tag expression
%    fn      filename with tags replaced by *



if nargin > 1
   tgr = parseParameter(varargin);
   texpr = tagExpressionToString(texpr, tagFromIndex(tgr, 1));
end

[~, tsplit, ~] = tagExpressionToTagNames(texpr);
fn = tsplit{1};
for i = 2:length(tsplit)
   fn = [fn, '*', tsplit{i}]; %#ok<AGROW>
end

end