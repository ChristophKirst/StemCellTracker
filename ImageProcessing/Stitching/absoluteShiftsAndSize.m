function [ashifts, asize] = absoluteShiftsAndSize(shifts, imagesizes)
%
% [ashifts, asize] = absoluteShiftsAndSize(shifts, imagesizes)
%
% description:
%     calculates the full size of the algined image and absolute shifts of
%     the individual images give the relative shifts w.r.t to the first image
%     and individual image sizes 

if ~iscell(imagesizes) || ~iscell(shifts) || any(size(imagesizes) ~= size(shifts))
   error('absoluteShiftsAndSize: inconsistent input');
end

asiz = size(shifts);
dim = length(shifts{1});

maxsize = zeros(1,dim);
minshifts = zeros(1,dim);

for i = 1:numel(shifts)
   for d = 1:dim
      if shifts{i}(d) + imagesizes{i}(d) >= maxsize(d)
         maxsize(d) = shifts{i}(d) + imagesizes{i}(d);
      end
   
      if shifts{i}(d) < minshifts(d)
         minshifts(d) = shifts{i}(d);
      end
   end
end

asize = maxsize - minshifts;

ashifts = repmat({zeros(1,dim)}, asiz);
for i = 1:numel(shifts)
   ashifts{i} = shifts{i} - minshifts;
end

end

