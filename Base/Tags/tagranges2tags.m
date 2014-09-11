function tags = tagranges2tags(tagranges)
%
% tags = tagranges2tags(tagranges)
%
% dewcription:
%   converts tagranges in to tags struct array 
%
% See also: tags2tagranges

names = fieldnames(tagranges);
nnames = length(names);

if isempty(names)
   tags = struct;
   return
end

% assemble tag values for each name
vals = cell(1, nnames);
for n = 1:nnames
   vals{n} = tagranges.(names{n});
end

grc = ndgridc(vals{:});

stc = cell(1, 2 * nnames);
for i = 1:nnames
   stc{2*i-1} = names{i};
   gg = grc{i};
   stc{2*i}   = gg(:);
end

tags = struct(stc{:});
   
end