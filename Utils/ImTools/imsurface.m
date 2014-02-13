function [vertices, faces, normals] = imsurface(labeledimage)
%
% [vertices, facesangles, normals] = imisosurface(labeledimage)
%
% description:
%    calculates surfaces (vertices, facesangles) and surface normals for objects
%    in a labeled image
%
% input:
%    labeledimage   labeled object's image
%
% output:
%    vertices         vertices of each object surface as cell array
%    faces            faces
%    normals          normals
%
% See also: imsurfaceplot3d, impixelsurface, isosurface, isonormals, patch

label = imlabel(labeledimage);
isize = size(labeledimage);
nlabel = length(label);

vertices = cell(nlabel);
faces = cell(nlabel);
if nargout == 3
   normals = cell(nlabel);
end

for i = 1:nlabel
   l = label(i);
   obj = labeledimage == l;
   
   % reduce calculation to bounding box
   [bmin, bmax] = imboundingbox(obj);
   bmin = max(bmin - 1, 1);
   bmax = min(bmax + 1, isize);
   obj = imextract(obj, bmin, bmax);
   
   % call isosurface / isonormals
   [f,v] = isosurface(obj, 0.5);
   if nargout == 3
      n = isonormals(obj, v);
   end
   
   % correct for x,y exchange and assign outputs
   v = v(:, [2 1 3]);  
   vertices{i} = v + repmat(bmin, size(v,1), 1) - 1;
   faces{i} = f;
   
   if nargout == 3
      normals{i} = n(:,[2 1 3]);
   end
   
end

end