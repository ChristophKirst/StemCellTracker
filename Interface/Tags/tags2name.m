function name = tags2name(tfrmt, varargin)
%
% name = tags2name(tfrmt, tags)
%
% description:
%    generates a name from a tagged format string tfrmt and tag values tags
%    replacements are done as follows:
%    for a field with name xxx in the struct tags the substring 
%       <xxx>    is replaced by the value, 
%       <xxx,k>  is replaced by the value xxx using k digits with trailling zeros
%       <xxx,s>  denotes a string 
%
% input:
%    tfrmt       the tagged string
%    tags        struct with the image specifications as tags.name -> val or parameter list: 'name', val, ...
%
% output:
%    name        the name with tag replaced by values
%
% See also: name2tags, tagformat, num2str0, tagformat2tagnames


tags = parseParameter(varargin{:});

if nargin < 2 || isempty(tags) || isemptystruct(tags)
   name = tfrmt;
   return 
end

[tnames, tsplit, tinfo] = tagformat2tagnames(tfrmt);

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
            res{tinfo(i).pos(k)} = num2str(tags.(tagnames{p}));
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
   
   



