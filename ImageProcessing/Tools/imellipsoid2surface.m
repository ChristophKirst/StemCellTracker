function [vertices, faces, normals] = imellipsoid2surface(centroids, mainaxes, param)
%
% [vertices, faces, normals] = imellipsoid2surface(centroids, mainaxes, param)
%
% description:
%    converts ellipsoids to surfaces
%
% input:
%    centroids   centrodis of ellipsoids
%    mainaxes    main axes of ellipsiods
%    param       (optional) parameter struct wirh entries
%                .points   number of vertices
%
% output:
%     vertices     vertices of each object surface as cell array
%     faces        faces
%     normals      (optional) normals

if nargin < 3
   param = [];
end

npoints = getParameter(param, 'points', 25);

if isempty(centroids)
   vertices = {};
   faces = {};
   normals = {};
   return
end

dim = size(centroids,1);
ne = size(centroids,2);

if dim ~= 3
   error('imellipsoid2surface: expects 3d elliposids')
end

if size(mainaxes,1) ~= dim*dim || ne ~= size(mainaxes,2)
   error('imellipsoid2surface: diemsnions mismtach')
end


% calcualte surfaces

vertices = cell(ne,1);
faces = cell(ne,1);
if nargout > 2
   normals = cell(ne,1);
end

npoints2 = (npoints + 1) * (npoints + 1);

for i = 1:ne
   ax = mainaxes(1:3,i);
   ay = mainaxes(4:6,i);
   az = mainaxes(7:9,i);
   sh = repmat(centroids(:,i), 1, npoints2);
   
   rx = sqrt(sum(ax.^2));
   ry = sqrt(sum(ay.^2));
   rz = sqrt(sum(az.^2));
   
   ax = ax / rx;
   ay = ay / ry;
   az = az / rz;
   
   [x,y,z] = ellipsoid(0, 0, 0, rx,ry,rz, npoints);
   faces{i} = convhull(x,y,z);
   
   if nargout > 2
      [nx,ny,nz] = surfnorm(x,y,z);
   end
   
   %rotate onto axes
   xr = kron(ax, x(:)');
   yr = kron(ay, y(:)');
   zr = kron(az, z(:)'); 
   
   vertices{i} = (xr + yr + zr + sh)';
   %vertices{i} = v(:,[1 2 3]);
   
   
   if nargout > 2
      nx = kron(ax, nx(:)');
      ny = kron(ay, ny(:)');
      nz = kron(az, nz(:)');   
      normals{i} = (nx + ny + nz + sh)';
      %normals{i} = n(:,[2 1 3]);
   end
end

end