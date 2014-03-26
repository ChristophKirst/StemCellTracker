function coords = immask2coords(mask)
%
% coords = immask2coords(mask)
%
% description:
%    returns coordinates of the non-zero pixel in mask
%
% input:
%    mask          bw mask
%
% output:
%    coords        pixel coordinates in the form [p, q (, l)] 
%
%
% See also: imcoords2mask, imsub2ind, poly2mask

coords = imind2sub(size(mask), find(mask));

end