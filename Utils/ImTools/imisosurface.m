function [xyz, tri, nrm] = imisosurface(labeledimage)
%
% [xyz, triangles, normals] = imisosurface(labeledimage)
%
% description:
%    calculates surfaces (xyz, triangles) and surface normals for objects
%    in a labeled image
%
% input:
%    labeledimage   labeled object's image
%
% output:
%    xyz            coordinates for each object as cell array
%    tri            triangulation faces of surface as cell array
%    nrm            normals
%
% See also: implot3dsurface, isosurface, isonormals, patch

label = imlabel(labeledimage);
isize = size(labeledimage);
nlabel = length(label);

xyz = cell(nlabel);
tri = cell(nlabel);
if nargout == 3
   nrm = cell(nlabel);
end

for i = 1:nlabel
   l = label(i);
   obj = labeledimage == l;
   [bmin, bmax] = imboundingbox(obj);
   bmin = max(bmin - 1, 1);
   bmax = min(bmax + 1, isize);
   obj = imextract(obj, bmin, bmax);
   [f,v] = isosurface(obj, 0.5);
   xyz{i} = v + repmat(bmin, size(v,1), 1) - 1;
   tri{i} = f;
   
   if nargout == 3
      nrm{i} = isonormals(obj, v);
   end
end

end