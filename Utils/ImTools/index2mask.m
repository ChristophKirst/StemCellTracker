function mask = index2mask(imagesize, indices)
%
% mask = index2mask(imagesize, points)
%
% description:
%    returns bw image with pixels at points set to white
%
% input:
%    imagesize     [h, w] of return mask
%    points        coordinates of pixels as rows
%
% output:
%    mask          bw mask
%
% See also: poly2mask

if length(imagesize) < 2
   imagesize = [imagesize imagesize];
end

mask = zeros(imagesize(1), imagesize(2));

mask(sub2ind(imagesize, indices(1,:), indices(2,:))) = 1;

end