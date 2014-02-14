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
% note: the triangulation is assumed tobe in matlab index format
%       and in pixel coordinates
%
% output:
%    surface        ISurface surface reference
%
% See also: imarisset

[imaris, varargin, nargin] = imarisvarargin(varargin);

if nargin < 3 || nargin > 5
   error('imarissetvolume: expect 3-6 input arguments');
end

add = 0;
if ischar(varargin{1})
   surfacename = varargin{1};
   surface = imarisgetobject(imaris, surfacename, 'Surface');
   if isempty(surface)
      surface = imariscreatesurface(imaris, surfacename);
      add = 1;
   end 
   varargin = varargin(2:end);
   nargin = length(varargin);
elseif isimaristype(imaris, varargin{1}, 'Surface')
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
   timepoint = varargin{5};
end

surface.RemoveAllSurfaces
for i = 1:nsurfaces
  %vSurfaceHull.AddSurface(xyz{i}(:,[2, 1, 3])-1, tri{i}-1,  nrm{i}(:,[2, 1, 3]), timepoint);
  
  % convert vertices to space coordinates
  vertices{i} = imarispixel2space(imaris, vertices{i});
  surface.AddSurface(vertices{i}, faces{i}-1, normals{i}, timepoint);
end

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