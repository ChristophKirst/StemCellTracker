function [surface, varargout] = impixelsurface(labeledimage, stats)
%
% surface = impixelsurface(labeledimage)
% [surface, stats] = impixelsurface(labeledimage, stats)
%
% description: 
%    returns the surface of the labeld image
%
% input:
%    label    the labeled image
%    stats    (optional) statistics of labeled image
%
% output:
%    surface  labeled pixels on the surface of the labeled objects
%    stats    (optional) updated statistics of labeld image
%
% See also: bwpixelsurface

if nargin < 2
   stats = [];
end

stats = imstatistics(labeledimage, stats, 'PixelIdxList', 'BoundingBox');

isize = size(labeledimage);
dim = length(isize);
if dim < 2 || dim > 3
   error('impixelsurface: expect 2d or 3d labeled image');
end

surface = zeros(isize);

for l = 1:length(stats)
   
   bbox = stats(l).BoundingBox;
   
   %trailling size 1 does not work with matlab !
   if bbox(dim)==bbox(2*dim)
       if bbox(2*dim) < isize(dim)
           bbox(2*dim) = bbox(dim) + 1;
       else
           bbox(dim) = bbox(dim) - 1;
           if bbox(dim) < 1
               bbox(dim) = 1;
           end
       end
   end
   
   bmin = bbox(1:dim);

   obj = imextract(labeledimage, bbox);
   obj = (obj ==l);
   objsurf = bwpixelsurface(obj, 'border');
   idxsurf = find(objsurf);
   
   if dim == 2
      
      [sx,sy] = ind2sub(size(obj), idxsurf);
      sx = sx + bmin(1) - 1;
      sy = sy + bmin(2) - 1;
      surface(sub2ind(isize, sx, sy)) = l;
      
   else
      
      [sx,sy,sz] = ind2sub(size(obj), idxsurf);
      sx = sx + bmin(1) - 1;
      sy = sy + bmin(2) - 1;
      sz = sz + bmin(3) - 1;
      surface(sub2ind(isize, sx, sy, sz)) = l;
      
   end
end

if nargout > 1
    varargout{1} = stats;
end

end