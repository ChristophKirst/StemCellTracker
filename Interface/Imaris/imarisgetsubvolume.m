function stack = imarisgetsubvolume(varargin)
%
% stack = imarisgetsubvolume(p0, q0, l0, dp, dq, dl, channel, timepoint)
% stack = imarisgetsubvolume([p0, q0, l0], [dp, dq, dl], channel, timepoint)
% stack = imarisgetsubvolume(object, ...)
% stack = imarisgetsubvolume(imaris, ...)
%
% description:
%   retrieves subset of the volumetric data from Imaris object at a given timpepoint and color channel
%
% input:
%    pql0           position in pixel coordinates
%    dpql           size
%    channel        channel number 
%    timepoint      timepoint number
%
%    object         Imarise IDataSet object
%    imaris         Imaris application instance
%
% output:
%    stack          volumetric data
%
% See also: imarisgetvolume

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

pql0 = varargin{1};
if length(pql0) == 1
   if nargin < 6 
      error('imarisgetsubvolume: expect 6 parameters: p0, q0, l0, dp, dq, dl');
   end
   pql0 = [varargin{1:3}];
   dpql = [varargin{1:3}];
   varargin = varargin(7:end);
   nargin = length(varargin);
   
elseif length(pql0) == 3
   dpql = [varargin{2}];
   if length(dpql) ~=3 
      error('imarisgetsubvolume: expect parameters: [p0, q0, l0], [dp, dq, dl]');
   end
   
   varargin = varargin(3:end);
   nargin = length(varargin);
   
else 
   error('imarissetsubvolume: input syntax error.');
   
end

% matlab first index = 1, Imaris first index = 0
pql0 = pql0 -1;
p0 = pql0(1); q0 = pql0(2); l0 = pql0(3);
dp = dpql(1); dq = dpql(2); dl = dpql(3);

if p0 < 0 || p0 > dataset.GetSizeX() - 1
   error('imarisget: starting position p0 out of bounds.');
end
p0 = uint32(p0);

if q0 < 0 || q0 > dataset.GetSizeY() - 1
   error('imarisgetsubvolume: starting position q0 out of bounds.');
end
q0 = uint32(q0);

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
stack = zeros([dp, dq, dl], datatype);

% Get the stack
switch char(dataset.GetType())
   case 'eTypeUInt8',
      % Java does not have unsigned ints
      arr = dataset.GetDataSubVolumeAs1DArrayBytes(p0, q0, l0, channel, timepoint, dp, dq, dl);
      stack(:) = typecast(arr, 'uint8');
   case 'eTypeUInt16',
      % Java does not have unsigned ints
      arr = dataset.GetDataSubVolumeAs1DArrayShorts(p0, q0, l0, channel, timepoint, dp, dq, dl);
      stack(:) = typecast(arr, 'uint16');
   case 'eTypeFloat',
      stack(:) = dataset.GetDataSubVolumeAs1DArrayFloats(p0, q0, l0, channel, timepoint, dp, dq, dl);
   otherwise,
      error('imarisgetsubvolume: Bad value for dataset.GetType().');
end

%data = permute(data,[2,1,3]);
end
