function ids = roiToImageIndex(ipos, isizes, rect)
%
% ids = roiToImageIndex(ipos, isizes, rect)
%
% description: find ids of the images that overlap with the rect given as [lowerleft, upperight] 
%
% input:
%      shifts     the image shifts
%      isizes     the image sizes
%      rect       the rectangle to check overlap with, ofthe from [p1, p2]  where p1 < p2 are row vectors of the corner coordinates

if isa(rect, 'ROI')
   rect = rect.boundingBox;
   rect = rect.toPixelArray();
end

n = numel(ipos);
ids = zeros(1, n) > 0;
for i = 1:n
   r = shiftsAndSizeToRect(ipos{i}, isizes{i});
   if isoverlapping(r, rect)
      ids(i) = 1;
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
   
   b = all(r1 <= r2);  
end