function imarisset(varargin)
%
% imarisset(stack, channel, timepoint)
% imarisset(vertices, faces, normals, timepoint)
% imarisset(imaris, ...)
% imarisset(objectname, ...)
% imarisset(imaris, objectname, ...)
%
% description:
%    set object in Imaris.
%
% input:
%    stack          image stack for volumetric data
%    vertices       vertices of surface triangulation
%    faces          faces of surface triagulation
%    normals        surface normals at vertices
%    channel        (optional) color channel (0)
%    timepoint      (optional) timepoint to put data (0)
%
%    imaris         Imaris application instance
%    objectname     name of Imaris object 
%
% See also: imarisget


function imarisput(varargin)
%
% imarisput(stack, channel, timepoint)
% imarisput(xyz, tri, nrm, timepoint)
%
% description:
%    sets a stack in Imaris in color channel and timepoint to stack
%    or setst a surface objects defined by vertices xyz, faces tri and
%    normals nrm (as cell arrays)
%
% input:
%    stack          image stack
%    channel        color channel to put data 
%    timepoint      timepoint to put data
%    xyz, tri, nrm  triangulation data with vertices xyz, faces tri and normals nrm
%
% See also: convexhulln

[imaris, varargin, nargin] = imarisvarargin(varargin{:});

if nargin < 1
   error('imarisput: wrong number of input arguments!')
end

if iscell(varargin{1})
   type = 'Surface';
   if nargin < 3
      error('imarisput: wrong number of input arguments!')
   end
   if nargin < 4
      varargin{4} = 0;
   end
else
   type = 'Volume';
   if nargin < 2 % volume
      varargin{2} = 0;
   end
   if length(varargin) < 3
      varargin{3} = 0;
   end
end


switch type
   case 'Volume'
      imarisputstack(imaris, varargin{1:3});
   case 'Surface'
      imarisputsurface(imaris, varargin{1:4});
   otherwise
      error('imarisput: output type unsupported!')
end

end




function imarisputsurface(vImarisApplication, xyz, tri, nrm, timepoint)
 
if ~iscell(xyz)
   xyz = {xyzX};
end

if ~iscell(tri)
   tri = {tri};
end
   
nSurfaces = length(xyz);
if nSurfaces ~= length(tri) 
   error('imarisput: surface sizes do not agree');
end
if nSurfaces ~= length(nrm) 
   error('imarisput: surface sizes do not agree');
end

vFactory = vImarisApplication.GetFactory;
vSurpassScene = vImarisApplication.GetSurpassScene;

vSurfaceHull = vFactory.CreateSurfaces;

for i = 1:nSurfaces

  %vSurfaceHull.AddSurface(xyz{i}(:,[2, 1, 3])-1, tri{i}-1,  nrm{i}(:,[2, 1, 3]), timepoint);
  vSurfaceHull.AddSurface(xyz{i}, tri{i}-1, nrm{i}, timepoint);
  
end

vSurfaceHull.SetName('Matlab Surface');
%vSurfaceHull.SetColorRGBA(vFilaments.GetColorRGBA);
vSurpassScene.AddChild(vSurfaceHull, -1);

end





[imaris, varargin, nargin] = imarisvarargin(varargin{:});

if nargin < 1
   error('imarisput: wrong number of input arguments!')
end

if iscell(varargin{1})
   type = 'Surface';
   if nargin < 3
      error('imarisput: wrong number of input arguments!')
   end
   if nargin < 4
      varargin{4} = 0;
   end
else
   type = 'Volume';
   if nargin < 2 % volume
      varargin{2} = 0;
   end
   if length(varargin) < 3
      varargin{3} = 0;
   end
end


switch type
   case 'Volume'
      imarisputstack(imaris, varargin{1:3});
   case 'Surface'
      imarisputsurface(imaris, varargin{1:4});
   otherwise
      error('imarisput: output type unsupported!')
end

end



function imarisputstack(mImarisApplication, stack, channel, timepoint)

fprintf('pushing stack');

% transform stack
%stack = permute(stack, [2 1 3]);

% create an alias
iDataSet = mImarisApplication.GetDataSet();

% check whether we have some voxels in the dataset
if isempty(iDataSet)
    
    % Create and store a new dataset
    sz = size(stack);
    if numel(sz) == 2
        sz = [sz 1];
    end
    iDataSet = createDataset(mImarisApplication, class(stack), sz(1), sz(2), sz(3), 1, 1);

end

% Convert channel and timepoint to 0-based indexing
%channel = channel - this.mIndexingStart;
%timepoint = timepoint - this.mIndexingStart;

% Check that the requested channel and timepoint exist
if channel > (iDataSet.GetSizeC() - 1)
    error('The requested channel index %g is out of bounds.', channel);
end
if timepoint > (iDataSet.GetSizeT() - 1)
    error('The requested time index %g is out of bounds.', timepoint);
end

% Get the dataset class
switch char(iDataSet.GetType())
    case 'eTypeUInt8',  outDatatype = 'uint8';
    case 'eTypeUInt16', outDatatype = 'uint16';
    case 'eTypeFloat',  outDatatype = 'single';
    otherwise,
        error('Bad value for iDataSet.GetType().');
end

% Check that the input and output datatypes match
if ~isa(stack, outDatatype)
    error('Data type mismatch. %s != %s', class(stack), outDatatype);
end

% Check that the size matches
outSizes = imarisgetsize(mImarisApplication);
if ismatrix(stack)
    sizes = [size(stack) 1];
else
    sizes = size(stack);
end
if any(sizes(1 : 3) ~= outSizes(1 : 3))
    disp(sizes);
    disp(outSizes);
    error('Data volume size mismatch.');
end

% Set the stack
switch char(iDataSet.GetType())
    case 'eTypeUInt8',   
        iDataSet.SetDataVolumeAs1DArrayBytes(stack(:), channel, timepoint);
    case 'eTypeUInt16',
        iDataSet.SetDataVolumeAs1DArrayShorts(stack(:), channel, timepoint);
    case 'eTypeFloat',
        iDataSet.SetDataVolumeAs1DArrayFloats(stack(:), channel, timepoint);
    otherwise,
        error('Bad value for iDataSet.GetType().');
end

end



function iDataset = createDataset(mImarisApplication, datatype, sizeX, sizeY, ...
     sizeZ, sizeC, sizeT, voxelSizeX, voxelSizeY, voxelSizeZ, deltaTime)
% IceImarisConnector:  createDataset (public method)
% 
% DESCRIPTION
% 
%   This method creates an Imaris dataset and replaces current one.
% 
% SYNOPSIS
%
%   (1) iDataset = createDataset(datatype, sizeX, sizeY, sizeZ, sizeC, sizeT)
%   (2) iDataset = createDataset(datatype, sizeX, sizeY, sizeZ, sizeC, sizeT, ...
%                                voxelSizeX, voxelsSizeY, voxelSizeZ, deltaTime)
% 
% INPUT
% 
%   datatype  : one of 'uint8', 'uint16', 'single', Imaris.tType.eTypeUInt8,
%               Imaris.tType.eTypeUInt16, Imaris.tType.eTypeFloat
%   sizeX     : dataset width
%   sizeY     : dataset height
%   sizeZ     : number of planes
%   sizeC     : number of channels
%   sizeT     : number of timepoints
%   voxelSizeX: (optional, default = 1) voxel size in X direction
%   voxelSizeY: (optional, default = 1) voxel size in Y direction
%   voxelSizeZ: (optional, default = 1) voxel size in Z direction
%   deltaTime : (optional, default = 1) time difference between consecutive
%               time points
% 
% OUTPUT
% 
%   iDataset  : created DataSet
%
% EXAMPLE
%
%    % Create a 2-by-3-by-2 stack
%    data(:, :, 1) = [ 11 12 13; 14 15 16 ];
%    data(:, :, 2) = [ 17 18 19; 20 21 22];
%    data = uint8(data);
%
%    % Create a dataset with sizeX = 3, sizeY = 2 and sizeZ = 2
%    conn.createDataset('uint8', 3, 2, 2, 1, 1);
%
%    % Copy data into the Imaris dataset
%    conn.setDataVolumeRM(data, 0, 0);
%
% REMARK
%
%   If you plan to set data volumes in column-major form, swap sizeX and
%   sizeY.

if nargin ~= 7 && nargin ~= 11
    % The this parameter is hidden
    error('6 or 10 input parameters expected.');
end

% Default voxel sizes
if nargin == 7
    voxelSizeX = 1;
    voxelSizeY = 1;
    voxelSizeZ = 1;
    deltaTime  = 1;
end


% Imaris datatype
switch char(datatype)
    case {'uint8', 'eTypeUInt8'},
        classDataSet = Imaris.tType.eTypeUInt8;
    case {'uint16', 'eTypeUInt16'},
        classDataSet = Imaris.tType.eTypeUInt16;
    case {'single', 'eTypeFloat'},
        classDataSet = Imaris.tType.eTypeFloat;
    otherwise,
        error('Bad data type.');
end

% Create the dataset
iDataset = mImarisApplication.GetFactory().CreateDataSet();
iDataset.Create(classDataSet, sizeX, sizeY, sizeZ, sizeC, sizeT);

% Apply the spatial calibration
iDataset.SetExtendMinX(0);
iDataset.SetExtendMinY(0);
iDataset.SetExtendMinZ(0);
iDataset.SetExtendMaxX(sizeX * voxelSizeX);
iDataset.SetExtendMaxY(sizeY * voxelSizeY);
iDataset.SetExtendMaxZ(sizeZ * voxelSizeZ);

% Apply the temporal calibration
iDataset.SetTimePointsDelta(deltaTime);

% Set the dataset in Imaris
mImarisApplication.SetDataSet(iDataset);

end




