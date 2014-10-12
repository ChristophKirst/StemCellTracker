function id = tag2ind(tagranges, tag)
%
% id = tag2ind(tagranges, tag)
%
% description:
%   converts tag to the index i using tagranges

si = tagrangesize(tagranges);

names = fieldnames(tagranges);
nnames = length(names);

nt = length(tag);

id = zeros(1,nt);

for k = 1:nt

   fac = 1;
   id(k) = 0;
   for i = 1:nnames
      tr = tagranges.(names{i});
      idi = find(ismember(tr, tag(k).(names{i})));
      
      id(k) = id(k) + fac * (idi-1);
      fac = fac * si(i);
   end
   
   id(k) = id(k) + 1;

end

end