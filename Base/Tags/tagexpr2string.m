function name = tagexpr2string(texpr, varargin)
%
% name = tagexpr2string(texpr, tags)
%
% description:
%    generates a name from a tagged expression string texpr and tag values tags
%    replacements are done as follows:
%    for a field with name xxx in the struct tags the substring 
%       <xxx>    is replaced by the value, 
%       <xxx,k>  is replaced by the value xxx using k digits with trailling zeros
%       <xxx,s>  denotes a string 
%
% input:
%    texpr       the tagged string
%    tags        struct with the image specifications as tags.name -> val or parameter list: 'name', val, ...
%
% output:
%    name        the name with tag replaced by values
%
% See also: name2tags, tagexpr, num2str0, tagexpr2tagnames

if nargin < 2
   name = texpr;
   return 
end

if nargin == 2
   tags = varargin{1};
else
   tags = parseParameter(varargin{:});
end
%tags


 if isempty(tags) || isemptystruct(tags)
   name = texpr;
   return 
 end 

lt = length(tags);
if lt > 1
   name = cell(1, lt);
   for i = 1:lt
      name{i} = tagexpr2string(texpr, tags(i));
   end
   name = name';
   return
end

% process single tag struct
[tnames, tsplit, tinfo] = tagexpr2tagnames(texpr);

tagnames = fieldnames(tags);

res = repmat({''}, length(tsplit)-1, 1);
for i = 1:length(tnames)
   p = find(ismember(tagnames, tnames{i}));
   if isempty(p)
      for k = 1:length(tinfo(i).tag)
         res{tinfo(i).pos(k)} = tinfo(i).tag{k};
      end
   else
      p = p(1);
      for k = 1:length(tinfo(i).tag)
         tw = tinfo(i).width(k);
         if tw == 0
            res{tinfo(i).pos(k)} = var2char(tags.(tagnames{p}));
         else
            res{tinfo(i).pos(k)} = num2str0(tags.(tagnames{p}), tw);
         end
      end
   end
end

name = tsplit{1};
for i = 1:length(res);
   name = [name,  res{i}, tsplit{i+1}]; %#ok<AGROW>
end


end
   
   



