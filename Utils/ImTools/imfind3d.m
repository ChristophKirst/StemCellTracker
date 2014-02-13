function [h,w,l] = imfind3d(image, coords)
%
% [h,w,l] = imfind3d(image)
%
% descrition:
%     finds the coordinates of nonzero values in the image
%
% input:
%     image      image
%     coords     (optional) 'hwl' pixel coordinates, 'xyz' space coordinates
%
% output:
%     h,w,l       pixel coordinates

idx = find(image);
[h,w,l] = ind2sub(size(image), idx);

if nargin > 1
   if ischar(coords) && strcmp(coords, 'xyz') % swap h and w
      ww = h; h = w; w = ww;
   end
end

end