function surface = impixelsurface(labeledimage)
%
% surface = impixelsurface(labeledimage)
%
% description: 
%    returns the surface of the labeld image
%
% input:
%    label    the labeled image
%
% output:
%    surface  labeled pixels on the surface of the labeled objects
%

label = imlabel(labeledimage);
isize = size(labeledimage);
if length(isize) ~= 3
   error('imsurface: expect 3d bw image')
end

surface = labeledimage;

for l = label
   obj = labeledimage == l;
   [bmin, bmax] = imboundingbox(obj);
   obj = imextract(obj, bmin, bmax);
   objsurf = bwsurface(obj);
   idxsurf = find(obj - objsurf);
   [sx,sy,sz] = ind2sub(size(obj), idxsurf);
   sx = sx + bmin(1) - 1;
   sy = sy + bmin(2) - 1;
   sz = sz + bmin(3) - 1;
   surface(sub2ind(isize, sx, sy, sz)) = 0;
end

end