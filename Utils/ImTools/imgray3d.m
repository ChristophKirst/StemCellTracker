function image = imgray3d(stack, normalize)
%
% image = imgray3D(stack, normalize)
%
% description: 
%     takes a stack of gray images and moves them into hxwx1xz format 
%     for functions such as montage etc
%
% input:
%     stack      h x w x z stack of h x w images
%     normalize  optional global noramalization (1 = noramlize, 0 = ignorre) 
%
% output:
%     ie       cropped image
%
% See also: montage

if nargin < 2
   normalize = 1;
end
if ndims(stack) ~= 3
   error('imgray3D: expect stack of grayscale images')
end

[w, h, z] = size(stack);
image = reshape(stack, [w h 1 z]);

if normalize
   image = image / max(image(:));
end
   