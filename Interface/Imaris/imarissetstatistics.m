function istat = imarissetstatistics(varargin)
%
% istat = imarissetstatistics(name, values)
% istat = imarissetsurface(objectname, ...)
% istat = imarissetsurface(object, ...)
% istat = imarissetsurface(imaris, ...)
%
% description:
%    set statistic values in Imaris
%
% input:
% 
% output:
%    istat   statistics
%
% See also: imarisset

[imaris, varargin, nargin] = imarisvarargin(varargin);

%%

obj = imarisgetcurrentobject();
istat = obj.GetStatistics();

nstats = obj.GetNumberOfSurfaces();



%%

data = ones(1,nstats);

vSize = nstats;

vfac = {'Surface', '', '', '1'}';
vFactorNames = {'Category', 'Channel', 'Collection', 'Time'}';

aSimilarityNames   = cell(vSize, 1);
aSimilarityUnits   = cell(vSize, 1);
aSimilarityFactors = cell(size(vfac, 1), vSize);
aSimilarityIds     = 1:vSize;
for j = 1:vSize
    aSimilarityNames{j}      = ['Test Statistics'];
    aSimilarityUnits{j}      = 'um';
    aSimilarityFactors(:, j) = vfac;
    aSimilarityIds(j)        = j;
end

aSimilarityValues  = data;
 

%%
obj.AddStatistics(aSimilarityNames, aSimilarityValues, ...
      aSimilarityUnits, aSimilarityFactors, vFactorNames, aSimilarityIds);


%%

try
obj.AddStatistics(aSimilarityNames, aSimilarityValues, ...
      aSimilarityUnits, {}, {}, aSimilarityIds);
catch
    disp oops
end
  
  
  %%








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
   timepoint = varargin{5};
end

psize = imarisgetsize(imaris);
extend = imarisgetextend(imaris);
fac = (extend(2,:) - extend(1,:)) ./ psize;


surface.RemoveAllSurfaces
for i = 1:nsurfaces
  %vSurfaceHull.AddSurface(xyz{i}(:,[2, 1, 3])-1, tri{i}-1,  nrm{i}(:,[2, 1, 3]), timepoint);
  
  % convert vertices to space coordinates
  % vertices{i} = imarispixel2space(imaris, vertices{i});
  verts = impixel2space(psize, extend, vertices{i});
  
  nrmls = normals{i};
  nrmls = nrmls .* repmat(fac, size(nrmls,1),1);

  surface.AddSurface(verts, faces{i}-1, nrmls , timepoint);
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