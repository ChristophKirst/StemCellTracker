function si = tagrangesize(tagrs)        
%
% si = tagrangesize(tagrs)        
%
% description:
%      calculates sizes of the tag ranges
%
% input:
%      tagrs   struct array of tag ragnes
%
% output:
%      si      sizes
%
% See also: tagexpr


if isempty(tagrs)
   si = 0;
   return
end

names = fieldnames(tagrs);
nnames = length(names);

si= zeros(1,nnames);

for i = 1:nnames
   si(i) = length(tagrs.(names{i}));
end

end