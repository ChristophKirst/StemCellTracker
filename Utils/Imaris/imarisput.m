function imarisput(varargin)
%
% imarisput(stack, channel, timepoint)
% imarisput(X, K, timepoint)
%
% description:
%    sets a stack in Imaris in color channel and timepoint to stack
%    or setst a surface object defined by coordinates X and facets K
%
% See also: convexhulln

if length(varargin) < 2
   error('imarisput: wrong number of input arguments!')
end
if length(varargin) < 3
   varargin{3} = 0;
end

[imaris, varargin] = imarisvarargin(varargin);

if isscalar(varargin{2})
   imarisputstack(imaris, varargin{1:3});
else
   imarisputsurface(imaris, varargin{1:3});
end

end



function imarisputstack(mImarisApplication, stack, channel, timepoint)

% Create an alias
iDataSet = mImarisApplication.GetDataSet();

% Check whether we have some voxels in the dataset
if isempty(iDataSet)
    
    % Create and store a new dataset
    sz = size(stack);
    if numel(sz) == 2
        sz = [sz 1];
    end
    iDataSet = createDataset(mImarisApplication, class(stack), sz(1), sz(2), sz(3), 1, 1);

end

% Convert channel and timepoint to 0-based indexing
channel = channel - this.mIndexingStart;
timepoint = timepoint - this.mIndexingStart;

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
    error('Data type mismatch.');
end

% Check that the size matches
outSizes = this.getSizes();
if ismatrix(stack)
    sizes = [size(stack) 1];
else
    sizes = size(stack);
end
if any(sizes(1 : 3) ~= outSizes(1 : 3))
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







function imarisputsurface(mImarisApplication, X, K, timepoint)
 
if ~iscell(X)
   X = {X};
end

if ~iscell(K)
   K = {K};
end
   
nSurfaces = length(X);
if nSurfaces ~= length(K) 
   error('imarisput: surface sizes do not agree');
end


vFactory = mImarisApplication.GetFactory;
vSurpassScene = vImarisApplication.GetSurpassScene;

vSurfaceHull = vFactory.CreateSurfaces;

for i = 1:nSurfaces

  vXYZ = X{i};
  vConvexHull = K{i};
  vNumberOfPoints = size(vXYZ, 1);
  
  vPoints = false(vNumberOfPoints, 1);
  vPoints(vConvexHull(:)) = true;
  vPoints = find(vPoints);
  vVertices = vXYZ(vPoints, :);

  % remap vertex indices to our selection
  % and reorder triangle vertices (clockwise to counter)
  vPointsMap = zeros(vNumberOfPoints, 1);
  vPointsMap(vPoints) = 1:numel(vPoints);
  vTriangles = vPointsMap(vConvexHull(:, [1, 3, 2])) - 1;

  % calculate normals (do not normalize them, imaris will do it)
  % follow rays from center to vertices
  vMean = mean(vVertices, 1);
  vNormals = [vVertices(:, 1) - vMean(1), vVertices(:, 2) - vMean(2), ...
    vVertices(:, 3) - vMean(3)];

  vSurfaceHull.AddSurface(vVertices, vTriangles, vNormals, timepoint);
  
end

vSurfaceHull.SetName('Matlab Surface');
%vSurfaceHull.SetColorRGBA(vFilaments.GetColorRGBA);
vSurpassScene.AddChild(vSurfaceHull, -1);

end
