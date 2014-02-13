function imsurfaceplot3d(vertices, faces, normals)
%
% h = imsurfaceplot3d(vertices, faces, normals)
% h = imsurfaceplot3d(labeledimage)
%
% description:
%    plots the surfaces obtained with imsurface
%    all coordinates are assumed to be in pixel coordinates hwl
%
% input:
%    vertices       coordinates for each object as cell array
%    faces          facesangulation faces of surface as cell array
%    normals        normals
%    labeledimage   surfaces are infered form the labeled image
%
% See also: imsurface


if nargin == 1
   isize = size(vertices);
   [vertices, faces, normals] = imsurface(vertices);
else
   isize = [];
end

if ~iscell(vertices)
   vertices = {vertices};
end
if ~iscell(faces)
   faces = {faces};
end

nlabel = length(vertices);
if nlabel ~= length(faces)
   error('imsurfaceplot3d: inconsistent input sizes');
end

if nargin < 3
   normals = cell(1,nlabel);
else
   if ~iscell(normals)
      normals = {normals};
   end
   if nlabel ~= length(normals)
      error('imsurfaceplot3d: inconsistent input sizes');
   end   
end


%plot the surfaces 

hold on
cm = colormap;
ncm = length(cm);

for i = 1:nlabel
   
   fv.vertices = vertices{i};
   fv.faces = faces{i};
   col = cm(round((i-1)/nlabel * (ncm-1))+1, :); 
   if isempty(normals{i})
      patch(fv, 'FaceColor', col ,'EdgeColor','none');
   else
      patch(fv, 'FaceColor', col ,'EdgeColor','none', 'VertexNormals', normals{i});
   end

end

% some standard nice settings
daspect([1 1 1]); view(3); axis tight

if ~isempty(isize)
   xlim([0, isize(1)]); ylim([0, isize(2)]); zlim([0, isize(3)])
end

camlight
xlabel('h'); ylabel('w'); zlabel('l');

end



