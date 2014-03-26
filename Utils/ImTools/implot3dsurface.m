function implot3dsurface(vertices, faces, normals, param)
%
% implot3dsurface(vertices, faces, normals, param)
% implot3dsurface(labeledimage, param)
%
% description:
%    plots the surfaces obtained with imsurface
%    all coordinates are assumed to be in pixel coordinates pql
%
% input:
%    vertices       coordinates for each object as cell array
%    faces          triangulation faces of surface as cell array
%    normals        normals
%    labeledimage   surfaces are inferred form the labeled image
%    param          (optional)    parameter struct with entries
%                   .boundary     'pql' to close surfaces in specified directions ('pql')
%                   .color.data   color according to this data ([] = random)
%                   .color.map    use this colormap (colormap)
%                   .color.scale  scale color data
%
% See also: imsurface, imlabelcolormap


if (nargin ==2 && isstruct(faces)) 
   param = faces;
elseif  (nargin ==3 && isstruct(normals)) 
   param = normals;
elseif nargin < 4
   param = [];
end

bd = getParameter(param, {'boundary'}, 'pql');

if nargin ==1 || (nargin ==2 && isstruct(faces))
   isize = size(vertices);
   [vertices, faces, normals] = imsurface(vertices, bd);
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

if ~exist('normals', 'var')
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
cm = imlabelcolormap(nlabel, param);

for i = 1:nlabel
   fv.vertices = vertices{i};
   fv.faces = faces{i};
   col = cm(i, :); 
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
xlabel('p'); ylabel('q'); zlabel('l');

end



