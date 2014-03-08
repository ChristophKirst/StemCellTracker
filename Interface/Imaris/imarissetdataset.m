function dataset = imarissetdataset(varargin)
%
% dataset = imarissetdataset(datatype, sizeX, sizeY, sizeZ, sizeC, sizeT)
% dataset = imarissetdataset(datatype, sizeX, sizeY, sizeZ, sizeC, sizeT, 
%                               voxelSizeX, voxelSizeY, voxelSizeZ, deltaTime)
% 
% description: 
%   create an Imaris IDataSet object and replace current one.
% 
% input:
%   datatype    'uint8', 'uint16', 'single', 
%               Imaris.tType.eTypeUInt8, Imaris.tType.eTypeUInt16, Imaris.tType.eTypeFloat
%   sizeX       dataset width
%   sizeY       dataset height
%   sizeZ       number of planes
%   sizeC       number of channels
%   sizeT       number of timepoints
%   voxelSizeX  (optional, default = 1) voxel size in X direction
%   voxelSizeY  (optional, default = 1) voxel size in Y direction
%   voxelSizeZ  (optional, default = 1) voxel size in Z direction
%   deltaTime   (optional, default = 1) time difference between consecutive time points
% 
% output:
%   dataset     new Imaris IDataSet
%
% See also: imarisset

[imaris, varargin, nargin] = imarisvarargin(varargin);

if nargin ~= 6 && nargin ~= 10
    error('imariscreatedataset: 6 or 10 input parameters expected.');
end

datatype = varargin{1};
sizeX    = varargin{2};
sizeY    = varargin{3};
sizeZ    = varargin{4}; 
sizeC    = varargin{5};
sizeT    = varargin{6};

if nargin == 6
    voxelSizeX = 1;
    voxelSizeY = 1;
    voxelSizeZ = 1;
    deltaTime  = 1;
else
    voxelSizeX = varargin{7};
    voxelSizeY = varargin{8};
    voxelSizeZ = varargin{9};
    deltaTime  = varargin{10}; 
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
        error('imariscreatedataset: bad data type %s.', type);
end

% Create the dataset
dataset = imaris.GetFactory().CreateDataSet();
dataset.Create(classDataSet, sizeX, sizeY, sizeZ, sizeC, sizeT);

% Apply the spatial calibration
dataset.SetExtendMinX(0);
dataset.SetExtendMinY(0);
dataset.SetExtendMinZ(0);
dataset.SetExtendMaxX(sizeX * voxelSizeX);
dataset.SetExtendMaxY(sizeY * voxelSizeY);
dataset.SetExtendMaxZ(sizeZ * voxelSizeZ);

% Apply the temporal calibration
dataset.SetTimePointsDelta(deltaTime);

% Set the dataset in Imaris
imaris.SetDataSet(dataset);

end
