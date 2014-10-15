function [tags, names, ttypes] = tagExpressionToTags(texpr, names, varargin)
%
% tags = tagExpressionToTags(texpr, name, param)
%
% description:
%    finds tags according to the tag expression texpr in name
%
% input:
%    tformat     tag expression string
%    names       string or cellstr, or [] = filenames that match texpr to infer tags from
%    param       (optional) paramter struct with entries
%                .tagnames    list with tag names to return in tags ([] = all)
%
% See also: tagexpr, tagexpr2tagnames, num2str0

% convert tagexpr to regexpr

[tnames, ~, ti] = tagExpressionToTagNames(texpr);
re = tagExpressionToRegularExpression(texpr);

param = parseParameter(varargin);
tagnames = getParameter(param, 'tagnames', []);

if isempty(tagnames)
   tagnames = tnames;
else
   if ~iscellstr(tagnames)
      if ischar(tagnames)
         tagnames = {tagnames};
      else
         error('tagExpressionToTags: expects cell or strings as third argument');
      end
   end
end
n = length(tnames);

%check = getParameter(param, 'check', false);
if nargin < 2 || isempty(names)
   names = tagExpressionToFiles(texpr);
   %check = false; % no need for checking as file list is in accordance with texpr
end

if ~iscell(names) 
   names = {names};
end

% multiple occurences of a tag
% if check && all(cellfun(@length, {ti.pos}) == 1)
%    check = false;
% end

%find tags
m = length(names);
%re
for i = m:-1:1
   %s = regexp(names{i}, re, 'names')
   tags(i) = regexp(names{i}, re, 'names');
end

% filter out required tagsnames
for i = 1:length(tagnames);
   k = find(~cellfun(@isempty, strfind(tnames,tagnames{i})));
   if isempty(k)
      error('tagExpressionToTags: unable to locate tag %s', tagnames{i});
   end
   ttypes{i} = ti(k(1)).type; %#ok<AGROW>
end
tags = rmfield(tags, setdiff(fieldnames(tags), tagnames));

% convert to types
if numel(tags) > 0
   for i = n:-1:1
      if isfield(tags, tagnames{i})
         if strcmp(ttypes{i}, 'd')
            dat = num2cell(str2double({tags.(tagnames{i})}));
            [tags.(tagnames{i})] = dat{:};
         end
      else
         error('tagExpressionToTags: unable to locate tag %s', tnames{i});
      end
   end
end

% %check for consistency if multiple occurences of a tag
% if check
%    names2 = tagExpressionToString(texpr, tags);
%    if ischar(names2)
%       names2 = {names2};
%    end
%   
% %    var2char(names)
% %    var2char(names2)
% %    class(names)
% %    class(names2)
%   
%    ids = ~cellfun(@isequal, names(:), names2(:));
%    
%    if total(ids) > 0
%       warning('tagExpressionToTags: some strings %s do not match tagexpr: %s!', var2char(names(ids)), texpr)
%       names = names(~ids);
%       tags  = tags(~ids);
%    end
% end










end
   



