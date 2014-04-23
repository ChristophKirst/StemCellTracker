function [tags, ttypes] = name2tags(tfrmt, name, tagnames)
%
% tags = name2tags(tfrmt, tagnames)
%
% description:
%    finds tags according to the tag format tfrmt in name
%
% input:
%    tformat     tag format string
%    name        name string or cell of strings
%    tagnames    (optional) list and ordering of tag names to be returned
%
% See also: tags2name, tagformat, tagformat2tagnames, num2str0

% convert tagformat to regexpr

[tnames, ~, ttype] = tagformat2tagnames(tfrmt);
re = tagformat2regexp(tfrmt);

if nargin == 3
   if ~iscellstr(tagnames)
      if ischar(tagnames)
         tagnames = {tagnames};
      else
         error('name2tags: expects cell or strings as third argument');
      end
   end
else
   tagnames = tnames;
end


n = length(tnames);
if ~iscell(name) 
   name = {name};
end

m = length(name);
for i = m:-1:1
   regexp(name{i}, re, 'match')
   tags(i) = regexp(name{i}, re, 'names');
end

for i = 1:length(tagnames);
   k = find(~cellfun(@isempty, strfind(tnames,tagnames{i})));
   if isempty(k)
      error('name2tags: unable to locate tag %s', tagnames{i});
   end
   ttypen{i} = ttype{k}; %#ok<AGROW>
end

tags = rmfield(tags, setdiff(fieldnames(tags), tagnames));

for i = n:-1:1
   if isfield(tags, tagnames{i})
      if strcmp(ttypen{i}, 'd')
         dat = num2cell(str2double({tags.(tagnames{i})}));
         [tags.(tagnames{i})] = dat{:};
      end
   else
      error('name2tags: unable to locate tag %s', tnames{i});
   end
end

if nargout > 1
   ttypes = ttypen;
end

end
   



