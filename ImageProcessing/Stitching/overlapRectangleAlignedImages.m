function ids = overlapRectangleAlignedImages(shifts, isizes, rect)
%
%
% description: find ids of the images that overlap with the rect given as [lowerleft, upperight] 

n = numel(shifts);
ids = zeros(1, n);
for i = 1:n
   r = shiftsAndSizeToRect(shifts{i}, isizes{i});
   if isoverlapping(r, rect)
      ids(n) = 1;
   end
end

end



function r = shiftsAndSizeToRect(shift, isize)
   r = [shift(:), shift(:) + isize(:) - 1];
end

function b = isoverlapping(a, b)
   a1 = a(:,1); a2 = a(:,2);
   b1 = b(:,1); b2 = b(:,2);
   
   r1 = max(a1, b1);
   r2 = min(a2, b2);
   
   b = ~any(r1 < r2);  
end