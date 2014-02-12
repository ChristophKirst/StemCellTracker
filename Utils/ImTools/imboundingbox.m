function [minpos, maxpos] = imboundingbox(image)
%
% [minpos, maxpos] = imboundingbox(image)
%
% description:
%    finds the boundin box in pixel coordinates of the bw image
%
% input:
%    image  bw image
%
% output:
%    minpos  lower corner of bbox
%    maxpos  upper corner of bbox


d = ndims(image);
indx =  find(image > 0);

if d == 2
   [ix, iy] = ind2sub(size(image), indx);
   ii = [ix iy];
else
   [ix, iy, iz] = ind2sub(size(image), indx);
   ii = [ix iy iz];
end

minpos = min(ii);
maxpos = max(ii);

end