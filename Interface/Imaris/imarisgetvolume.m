function stack = imarisgetvolume(varargin)
%
% stack = imarisgetvolume(channel, timepoint)
% stack = imarisgetvolume(object, ...)
% stack = imarisgetvolume(imaris, ...)
%
% description:
%   retrieves volumetric data from Imaris object at a given timpepoint and color channel
%
% input:
%    channel  : channel number 
%    timepoint: timepoint number
%
%    object         Imarise IDataSet object
%    imaris         Imaris application instance
%
% output:
%    stack          volumetric data
%
% note: The IVolume object in Imaris is just a wrapper, all data is stored in DataSet
%
% See also: imarisset

[imaris, varargin, nargin] = imarisvarargin(varargin);

if nargin > 3
   error('imarissetvolume: expect 0-3 input arguments');
end

if isimaristype(imaris, varargin{1}, 'DataSet')
   dataset = varargin{1};
   varargin = varargin(2:end);
   nargin = length(varargin);
else
   dataset = imarisgetdataset(imaris);
end

% Check whether we have some voxels in the dataset
if isempty(dataset) || dataset.GetSizeX() == 0
   stack = [];
   return
end


if nargin < 1
   channel = 0;
else
   channel = varargin{1};
end
%channel = channel - 1;
if channel > (dataset.GetSizeC() - 1)
   error('imarissetvolume: color index is out of bounds.');
end

if nargin < 2
   timepoint = 0;
else
   timepoint = varargin{2};
end
%timepoint = timepoint - 1;
if timepoint > (dataset.GetSizeT() - 1)
   error('imarissetvolume: time index is out of bounds.');
end


% Get the dataset class
switch char(dataset.GetType())
   case 'eTypeUInt8',   datatype = 'uint8';
   case 'eTypeUInt16',  datatype = 'uint16';
   case 'eTypeFloat',   datatype = 'single';
   otherwise,
      error('imarissetvolume: Bad value for dataset.GetType().');
end

% Allocate memory
stack = zeros([dataset.GetSizeX(), dataset.GetSizeY(), ...
   dataset.GetSizeZ()], datatype);

% Get the stack
switch char(dataset.GetType())
   case 'eTypeUInt8',
      % Java does not have unsigned ints
      arr = dataset.GetDataVolumeAs1DArrayBytes(channel, timepoint);
      stack(:) = typecast(arr, 'uint8');
   case 'eTypeUInt16',
      % Java does not have unsigned ints
      arr = dataset.GetDataVolumeAs1DArrayShorts(channel, timepoint);
      stack(:) = typecast(arr, 'uint16');
   case 'eTypeFloat',
      stack(:) = dataset.GetDataVolumeAs1DArrayFloats(channel, timepoint);
   otherwise,
      error('imarissetvolume: Bad value for dataset.GetType().');
end

%stack = permute(stack,[2,1,3]);

end