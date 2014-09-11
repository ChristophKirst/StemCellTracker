function [tags, names, ttypes] = tagexpr2tags(texpr, names, varargin)
%
% tags = tagexpr2tags(texpr, name, param)
%
% description:
%    finds tags according to the tag expression texpr in name
%
% input:
%    tformat     tag expression string
%    names       string or cellstr, or [] = filenames that match texpr to infer tags from
%    param       (optional) paramter struct with entries
%                .tagnames    list with tag names to return in tags ([] = all)
%                .check       weather to check for consistency
%
%
% See also: tagexpr, tagexpr2tagnames, num2str0

% convert tagexpr to regexpr

[tnames, ~, ti] = tagexpr2tagnames(texpr);
re = tagexpr2regexp(texpr);

param = parseParameter(varargin);
tagnames = getParameter(param, 'tagnames', []);

if isempty(tagnames)
   tagnames = tnames;
else
   if ~iscellstr(tagnames)
      if ischar(tagnames)
         tagnames = {tagnames};
      else
         error('tagexpr2tags: expects cell or strings as third argument');
      end
   end
end
n = length(tnames);

check = getParameter(param, 'check', true);
if nargin < 2 || isempty(names)
   names = tagexpr2files(texpr);
   check = false; % no need for checking as file list is in accordance with texpr
end

if ~iscell(names) 
   names = {names};
end

% multiple occurences of a tag
if check && all(cellfun(@length, {ti.pos}) == 1)
   check = false;
end

%find tags
m = length(names);
for i = m:-1:1
   %regexp(name{i}, re, 'name')
   tags(i) = regexp(names{i}, re, 'names');
end

%check for consistency if multiple occurences of a tag
if check
   names2 = tagexpr2string(texpr, tags);
   if length(tags) ==1
      names2 = {names2};
   end
   
   names
   names2
  
   ids = ~cellfun(@isequal, names(:), names2(:));
   
   if total(ids) > 0
      warning('tagexpr2tags: some strings %s do not match tagexpr: %s!', var2char(names(ids)), texpr)
      names = names(~ids);
      tags  = tags(~ids);
   end
end



% filter out required tagsnames
for i = 1:length(tagnames);
   k = find(~cellfun(@isempty, strfind(tnames,tagnames{i})));
   if isempty(k)
      error('tagexpr2tags: unable to locate tag %s', tagnames{i});
   end
   ttypes{i} = ti(k(1)).type; %#ok<AGROW>
end
tags = rmfield(tags, setdiff(fieldnames(tags), tagnames));

if numel(tags) > 0
   % convert to types
   for i = n:-1:1
      if isfield(tags, tagnames{i})
         if strcmp(ttypes{i}, 'd')
            dat = num2cell(str2double({tags.(tagnames{i})}));
            [tags.(tagnames{i})] = dat{:};
         end
      else
         error('tagexpr2tags: unable to locate tag %s', tnames{i});
      end
   end
end



end
   



