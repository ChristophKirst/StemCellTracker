function implot3dellipsoid(centroids, mainaxes, param)
%
% implot3dellipsoid(centroids, mainaxes, param)
% implot3dellipsoid(labeledimage, param)
%
% description:
%    plots the ellipsoids obtained with imellipsoid
%    all coordinates are assumed to be in pixel coordinates pql
%
% input:
%    centroids      centroids of the ellipsoids as cell array or  matrix
%    mainaxes       facesangulation faces of surface as cell array
%    normals        normals
%    labeledimage   surfaces are inferred form the labeled image
%    param          (optional)    parameter struct with entries
%                   .n            number of surface points
%                   .color        color parameter struct as in imlabelcolormap
%
% See also: imsurface, imlabelcolormap

if nargin < 3
   param = [];
end

[vertices, faces, normals] = imellipsoid2surface(centroids, mainaxes, param);
implot3dsurface(vertices, faces, normals, param)

end