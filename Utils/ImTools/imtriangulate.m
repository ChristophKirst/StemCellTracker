function tri = imtriangulate(labeledimage)
%
% tri = imtriangulate(labeledimage)
%
%
% description:
%     triangulate the surface of the labeled objects
%
% input:
%     labeledimage   labled image
%
% output:
%     tri
%
% See also: imsurface, bwsurface

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
   
   

   
   
   
   
   
end

end


surface = imsurface(labeledimage);


