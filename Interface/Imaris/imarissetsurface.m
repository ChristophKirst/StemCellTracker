function surface = imarissetsurface(varargin)
%
% surface = imarissetsurface(vertices, faces, normals, timepoint)
% surface = imarissetsurface(objectname, ...)
% surface = imarissetsurface(object, ...)
% surface = imarissetsurface(imaris, ...)
%
% description:
%    set surface in Imaris.
%
% input:
%    vetices, faces, normals  surface triangulation data
%    timepoint      (optional) timepoint to put data (0)
%
%    objectname     name of Imaris surface
%    object         Imarise ISurface surface
%    imaris         Imaris application instance
%
% output:
%    surface        ISurface surface reference
%
% note: the triangulation is assumed to be in matlab index format (starting at 1)
%       and in pixel coordinates
%
% See also: imarissetvolume

[imaris, varargin, nargin] = imarisvarargin(varargin);

if nargin < 3 || nargin > 5
   error('imarissetvolume: expect 3-6 input arguments');
end

add = 0;
if ischar(varargin{1})
   surfacename = varargin{1};
   surface = imarisgetobject(imaris, surfacename, 'Surfaces');
   if isempty(surface)
      surface = imariscreatesurface(imaris, surfacename);
      add = 1;
   end 
   varargin = varargin(2:end);
   nargin = length(varargin);
elseif isimaristype(imaris, varargin{1}, 'Surfaces')
   surface = varargin{1};
   varargin = varargin(2:end);
   nargin = length(varargin);
else
   surface = imariscreatesurface(imaris, 'MSurface');
   add = 1;
end

if nargin < 3
   error('imarissetvolume: expect 3-6 input arguments');
end
vertices = varargin{1};
faces = varargin{2};
normals = varargin{3};

if ~iscell(vertices)
   vertices = {vertices};
end
if ~iscell(faces)
   faces = {faces};
end
if ~iscell(normals)
   normals = {normals};
end

nsurfaces = length(vertices);
if nsurfaces ~= length(faces) 
   error('imarisput: surface parameter sizes do not agree');
end
if nsurfaces ~= length(normals) 
   error('imarisput: surface parameter sizes do not agree');
end

if nargin < 4
   timepoint = 0;
else
   timepoint = varargin{4};
end

psize = imarisgetsize(imaris);
extend = imarisgetextend(imaris);
fac = (extend(2,:) - extend(1,:)) ./ psize;

surface.RemoveAllSurfaces

% join data to list
nverts = cellfun(@length, vertices);
nfaces = cellfun(@length, faces);
nsurfaces = length(vertices);

vertices = cell2mat(vertices);
faces    = cell2mat(faces) - 1;
normals  = cell2mat(normals);

vertices =  impixel2space(psize, extend, vertices);
normals = normals .* repmat(fac, size(normals,1),1);
timepoint = timepoint * ones(1, nsurfaces);

surface.AddSurfacesList(vertices, nverts, faces, nfaces, normals, timepoint);

% for i = 1:nsurfaces
%   %vSurfaceHull.AddSurface(xyz{i}(:,[2, 1, 3])-1, tri{i}-1,  nrm{i}(:,[2, 1, 3]), timepoint);
%   
%   % convert vertices to space coordinates
%   % vertices{i} = imarispixel2space(imaris, vertices{i});
%   verts = impixel2space(psize, extend, vertices{i});
%   
%   nrmls = normals{i};
%   nrmls = nrmls .* repmat(fac, size(nrmls,1),1);
% 
%   surface.AddSurface(verts, faces{i}-1, nrmls , timepoint);
%   %surface.AddSurface(verts, faces{i}-1, nrmls);
% end

%surface.SetColorRGBA(xxx);

if add
   scene = imaris.GetSurpassScene;
   scene.AddChild(surface, -1);
end

end


% helper
function surface = imariscreatesurface(imaris, surfacename)
   factory = imaris.GetFactory;
   surface = factory.CreateSurfaces;
   surface.SetName(surfacename);
end