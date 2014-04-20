function tags = name2tags(tfrmt, name, tagnames)
%
% tags = name2tags(tfrmt, tagnames)
%
% description:
%    finds tags according to the tag format tfrmt in name
%
% input:
%    tformat     tag format string
%    name        name sting
%    tagnames    (optional) list and ordering of tag names to be returned
%
% See also: tags2name, tagformat, tagformat2tagnames, num2str0

% convert tagformat to regexpr

[tnames, tagw] = tagformat2tagnames(tfrmt);

re = tfrmt;
for i = 1:length(tnames)
   if tagw(i) == 0 
      repl = regexp(tfrmt, ['<\s*?' tnames{i} '\s*?>'],'match');
      if length(repl) ~= 1
         error('tags2name: cannot infer tag values')
      end
      re = strrep(re, repl{1}, ['(?<' tnames{i} '>\d+?)']);
   else
      repl = regexp(tfrmt,  ['<\s*?' tnames{i}  '\s*?,\s*?' num2str(tagw(i)) '\s*?>'], 'match');
      if length(repl) ~= 1
         error('tags2name: cannot infer tag values')
      end
      re = strrep(re, repl{1}, ['(?<' tnames{i} '>\d{' num2str(tagw(i)) '})']);
   end
end

if nargin == 3
   if ~iscellstr(tagnames)
      if ischar(tagnames)
         tagnames = {tagnames};
      else
         error('name2tags: expects cell or strings as third argument');
      end
   end   
   tnames = tagnames;
end


n = length(tnames);
if ~iscell(name) 
   name = {name};
end

m = length(name);
for i = m:-1:1
   tagsre(i) = regexp(name{i}, re, 'names');
end

tags = zeros(n,m);
for i = n:-1:1
   if isfield(tagsre, tnames{i})
      tags(i,:) = str2double({tagsre.(tnames{i})});
   else
      error('name2tags: unable to locate tag %s', tnames{i});
   end
end

end
   



