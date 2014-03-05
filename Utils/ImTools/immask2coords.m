function coords = immask2coords(mask)
%
% coords = index2mask(mask)
%
% description:
%    returns bw image with pixels at points set to white
%
% input:
%    mask          bw mask
%
% output:
%    coords        pixel coordinates in the form [h, w (, l)] 
%
%
% See also: imcoords2mask, imsub2ind, poly2mask

coords = imind2sub(size(mask), find(mask));

end