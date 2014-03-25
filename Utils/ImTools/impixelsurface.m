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
   bmin = bbox(1:dim);
   obj = imextract(obj, bbox);
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

end