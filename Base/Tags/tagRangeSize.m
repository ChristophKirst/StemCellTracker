function si = tagRangeSize(tagRange)        
%
% si = tagRangeSize(tagrs)        
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


if isempty(tagRange)
   si = 0;
   return
end

names = fieldnames(tagRange);
nnames = length(names);

si= zeros(1,nnames);

for i = 1:nnames
   si(i) = length(tagRange.(names{i}));
end

end