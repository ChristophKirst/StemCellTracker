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
dim = length(isize);
if dim < 2 || dim > 3
   error('impixelsurface: expect 2d or 3d bw image');
end
surface = zeros(isize);

for l = label
   
   obj = labeledimage == l;
   [bmin, bmax] = imboundingbox(obj);
   obj = imextract(obj, bmin, bmax);
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