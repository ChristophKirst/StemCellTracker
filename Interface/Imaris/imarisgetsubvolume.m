function stack = imarisgetsubvolume(varargin)
%
% stack = imarisgetsubvolume(h0, w0, l0, dh, dw, dl, channel, timepoint)
% stack = imarisgetsubvolume([h0, w0, l0], [dh, dw, dl], channel, timepoint)
% stack = imarisgetsubvolume(object, ...)
% stack = imarisgetsubvolume(imaris, ...)
%
% description:
%   retrieves subset of the volumetric data from Imaris object at a given timpepoint and color channel
%
% input:
%    hwl0           position in pixel coordinates
%    dhwl           size
%    channel        channel number 
%    timepoint      timepoint number
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

if nargin < 2 || nargin > 9
   error('imarisgetsubvolume: expect 2-9 input arguments');
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


if nargin < 2 || nargin > 8
   error('imarisgetsubvolume: expect 2-8 input arguments');
end

hwl0 = varargin{1};
if length(hwl0) == 1
   if nargin < 6 
      error('imarisgetsubvolume: expect 6 parameters: h0, w0, l0, dh, dw, dl');
   end
   hwl0 = [varargin{1:3}];
   dhwl = [varargin{1:3}];
   varargin = varargin(7:end);
   nargin = length(varargin);
   
elseif length(hwl0) == 3
   dhwl = [varargin{2}];
   if length(dhwl) ~=3 
      error('imarisgetsubvolume: expect parameters: [h0, w0, l0], [dh, dw, dl]');
   end
   
   varargin = varargin(3:end);
   nargin = length(varargin);
   
else 
   error('imarissetsubvolume: input syntax error.');
   
end

% matlab first index = 1, Imaris first index = 0
hwl0 = hwl0 -1;
h0 = hwl0(1); w0 = hwl0(2); l0 = hwl0(3);
dh = dhwl(1); dw = dhwl(2); dl = dhwl(3);

if h0 < 0 || h0 > dataset.GetSizeX() - 1
   error('imarisget: starting position h0 out of bounds.');
end
h0 = uint32(h0);

if w0 < 0 || w0 > dataset.GetSizeY() - 1
   error('imarisgetsubvolume: starting position w0 out of bounds.');
end
w0 = uint32(w0);

if l0 < 0 || l0 > dataset.GetSizeZ() - 1
   error('imarisgetsubvolume: starting position l0 out of bounds.');
end
l0 = uint32(l0);

 
if nargin < 1
   channel = 0;
else
   channel = varargin{1};
end
%channel = channel - 1;
if channel > (dataset.GetSizeC() - 1)
   error('imarisgetsubvolume: color index is out of bounds.');
end

if nargin < 2
   timepoint = 0;
else
   timepoint = varargin{2};
end
%timepoint = timepoint - 1;
if timepoint > (dataset.GetSizeT() - 1)
   error('imarisgetsubvolume: time index is out of bounds.');
end

% Get the dataset class
switch char(dataset.GetType())
   case 'eTypeUInt8',   datatype = 'uint8';
   case 'eTypeUInt16',  datatype = 'uint16';
   case 'eTypeFloat',   datatype = 'single';
   otherwise,
      error('imarisgetsubvolume: Bad value for IDataSet::GetType().');
end

% Allocate memory
stack = zeros([dh, dw, dl], datatype);

% Get the stack
switch char(dataset.GetType())
   case 'eTypeUInt8',
      % Java does not have unsigned ints
      arr = dataset.GetDataSubVolumeAs1DArrayBytes(h0, w0, l0, channel, timepoint, dh, dw, dl);
      stack(:) = typecast(arr, 'uint8');
   case 'eTypeUInt16',
      % Java does not have unsigned ints
      arr = dataset.GetDataSubVolumeAs1DArrayShorts(h0, w0, l0, channel, timepoint, dh, dw, dl);
      stack(:) = typecast(arr, 'uint16');
   case 'eTypeFloat',
      stack(:) = dataset.GetDataSubVolumeAs1DArrayFloats(h0, w0, l0, channel, timepoint, dh, dw, dl);
   otherwise,
      error('imarisgetsubvolume: Bad value for dataset.GetType().');
end

%data = permute(data,[2,1,3]);
end
