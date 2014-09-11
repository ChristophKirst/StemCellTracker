function fn = tagexpr2filename(texpr)
%
%  fn = tagexpr2filename(texpr)
%
% description:
%    returns filename pattern that is in accordance with the tag format
%
% input:
%    texpr   tag expression
%    fn      filename with tags replaced by *

[~, tsplit, ~] = tagexpr2tagnames(texpr);
fn = tsplit{1};
for i = 2:length(tsplit)
   fn = [fn, '*', tsplit{i}]; %#ok<AGROW>
end

end