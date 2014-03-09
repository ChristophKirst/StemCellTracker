function mask = imcoords2mask(imagesize, coords)
%
% mask = index2mask(imagesize, coords)
%
% description:
%    returns bw image with pixels at points set to white
%
% input:
%    imagesize     image size of return mask
%    coords        pixel coordinates in the form [p, q (, l)] 
%
% output:
%    mask          bw mask
%
% See also: imsub2ind, poly2mask

mask = zeros(imagesize);
mask(imsub2ind(imagesize, coords)) = 1;

end