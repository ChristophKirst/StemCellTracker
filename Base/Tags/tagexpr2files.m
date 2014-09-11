function fns = tagexpr2files(texpr, varargin)
%
% fns = tagexpr2files(texpr, param)
%
% description:
%      returns filnames that are in accordance with the tag format
%
% input:
%      texpr   tag expression
%      param   (optional) parameter struct with entries
%              .check    check for consistent filenames if multi tags
%
% output:
%      fns     filnames that match texpr
%
% See also: tagexpr


param = parseParameter(varargin);

re = tagexpr2regexp(texpr);
fns = tagexpr2filename(texpr);

fns = dirr(fns);

% sort out matching regexp
fns = fns(cellfun(@(x) ~isempty(regexp(x, re, 'once')), fns, 'UniformOutput', true));

if getParameter(param, 'check', true);
   %take care of multiple occurences of a tag
   [~,~,tinfo] = tagexpr2tagnames(texpr);
   
   if any(cellfun(@length, {tinfo.pos}) > 1)
      tags = tagexpr2tags(texpr, fns, 'check', false);
      fns2 = tagexpr2string(texpr, tags);
      fns = fns(cellfun(@isequal, fns, fns2));
   end
end

end