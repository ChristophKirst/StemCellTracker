function implotellipsoid(centroids, mainaxes, param)
%
% implotellipsoid(centroids, mainaxes, param)
% implotellipsoid(labeledimage, param)
%
% description:
%    plots the ellipsoids obtained with imellipsoid in 2d or 3d
%    all coordinates are assumed to be in pixel coordinates pql
%
% input:
%    centroids      centroids of the ellipsoids as cell array or  matrix
%    mainaxes       facesangulation faces of surface as cell array
%    normals        normals
%    labeledimage   surfaces are inferred form the labeled image
%    param          (optional)    parameter struct with entries
%                   .points       number of surface points
%                   .color        color parameter struct as in imlabelcolormap
%
% See also: imellipsoid, imlabelcolormap, implot3dsurface

if isempty(centroids)
   return
end

if nargin < 3
   param = [];
end

if iscell(centroids)
   centroids = cell2mat(centroids);
end
if iscell(mainaxes)
   mainaxes = cell2mat(mainaxes);
end

dim = length(centroids(:, 1))

if length(mainaxes(:, 1)) ~= 2*dim
   error('implotellipsoid: ellipsoids have inconsistent dimensions!')
end

if dim < 2 || dim > 3
   error('implotellipsoid: ellipsoids data not 2d or 3d!')
end


if dim == 3
   implot3dellipsoid(centroids, mainaxes, param);
else
  
   hold on
   
   nellipses = length(centroids(1,:));
   cm = imlabelcolormap(nellipses, param);
   
   n = getParameter(param, 'points', 20);
   t = 0:(2 * pi / n):(2 * pi);

   for i = 1:nellipses
      phi = atan2(mainaxes(2,i), mainaxes(1,i));
      a = sqrt(sum(mainaxes(1:2,i).^2));
      b = sqrt(sum(mainaxes(3:4,i).^2)); 
      
      x = centroids(1,i)+a * cos(phi) * cos(t)-b * sin(phi) * sin(t);
      y = centroids(2,i)+a * sin(phi) * cos(t)+b * cos(phi) * sin(t);
      
      
      col = cm(i, :); 
      plot(x,y, 'Color', col)
   end
  
end





