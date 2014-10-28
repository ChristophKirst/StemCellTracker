function [ashifts, asize] = absoluteShiftsAndSize(shifts, isizes)
%
% [ashifts, asize] = absoluteShiftsAndSize(shifts, isizes)
%
% description:
%     calculates the full size of the algined image and absolute shifts of
%     the individual images such that all shifts are minimally non-negative
%     the output can be used to assemble the image
%
% input:
%    shifts    image shifts
%    isizes    image sizes
% 
% output:
%    ashifts    absolte image shifts, i.e. no negative shifts
%    asize      the image size needed to place allimages 
%    

%shifts
%isizes

if ~iscell(isizes) || ~iscell(shifts) || numel(isizes) ~= numel(shifts)
%    var2char(isizes)
%    var2char(shifts)
   error('absoluteShiftsAndSize: inconsistent input');
end

asiz = size(shifts);
dim = length(shifts{1});

maxsize = shifts{1};
minshifts = maxsize;

for i = 1:numel(shifts)
   for d = 1:dim
      if shifts{i}(d) + isizes{i}(d) >= maxsize(d)
         maxsize(d) = shifts{i}(d) + isizes{i}(d);
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

