function dataset = imarissetsubvolume(varargin)
%
% stack = imarissetsubvolume(stack, h0, w0, l0, channel, timepoint)
% stack = imarissetsubvolume(stack, [h0, w0, l0], channel, timepoint)
% stack = imarissetsubvolume(object, ...)
% stack = imarissetsubvolume(imaris, ...)
%
% description:
%    tries to replace a sub volume in the active IDataSet of Imaris
%
% input:
%    hwl0           position in pixel coordinates
%    stack          subset data
%    channel        channel number 
%    timepoint      timepoint number
%
%    object         Imarise IDataSet object
%    imaris         Imaris application instance
%
% See also: imarisset

[imaris, varargin, nargin] = imarisvarargin(varargin);

if nargin < 2 || nargin > 7
   error('imarissetsubvolume: expect 2-7 input arguments');
end

if isimaristype(imaris, varargin{1}, 'DataSet')
   dataset = varargin{1};
   varargin = varargin(2:end);
   nargin = length(varargin);
else
   dataset = imarisgetdataset(imaris);
end

if nargin < 2 || nargin > 6
   error('imarissetsubvolume: expect 2-6 input arguments');
end

stack = varargin{1};


hwl0 = varargin{2};
if length(hwl0) == 1
   if nargin < 4 
      error('imarissetsubvolume: expect 3 parameters: h0, w0, l0');
   end
   hwl0 = [varargin{2:4}];
   varargin = varargin(5:end);
   nargin = length(varargin);
   
elseif length(hwl0) == 3
   
   varargin = varargin(3:end);
   nargin = length(varargin);
   
else 
   error('imarissetsubvolume: input syntax error.');
   
end

% matlab first index = 1, Imaris first index = 0
hwl0 = hwl0 -1;
h0 = hwl0(1); w0 = hwl0(2); l0 = hwl0(3);

if h0 < 0 || h0 > dataset.GetSizeX() - 1
   error('imarisget: starting position h0 out of bounds.');
end
h0 = uint32(h0);

if w0 < 0 || w0 > dataset.GetSizeY() - 1
   error('imarissetsubvolume: starting position w0 out of bounds.');
end
w0 = uint32(w0);

if l0 < 0 || l0 > dataset.GetSizeZ() - 1
   error('imarissetsubvolume: starting position l0 out of bounds.');
end
l0 = uint32(l0);

 
[dh, dw, dl] = size(stack);

dsize = imarisgetsize(imaris, dataset);

if any(hwl0 + [dh,dw,dl] > dsize -1)
   error('imarissetsubvolume: subvolume to big.');
end



if nargin < 1
   channel = 0;
else
   channel = varargin{1};
end
%channel = channel - 1;
if channel > (dataset.GetSizeC() - 1)
   error('imarissetsubvolume: color index is out of bounds.');
end

if nargin < 2
   timepoint = 0;
else
   timepoint = varargin{2};
end
%timepoint = timepoint - 1;
if timepoint > (dataset.GetSizeT() - 1)
   error('imarissetsubvolume: time index is out of bounds.');
end


% Get the dataset class
switch char(dataset.GetType())
   case 'eTypeUInt8',   datatype = 'uint8';
   case 'eTypeUInt16',  datatype = 'uint16';
   case 'eTypeFloat',   datatype = 'single';
   otherwise,
      error('imarissetsubvolume: Bad value for IDataSet::GetType().');
end

% Check that the input and output datatypes match
if ~isa(stack, datatype)
    warning('imarissetsubvolume: data type mismatch. %s != %s', class(stack), datatype);
    stack = cast(stack, datatype);
end

% Set the stack
switch char(dataset.GetType())
    case 'eTypeUInt8',   
        dataset.SetDataSubVolumeAs1DArrayBytes(stack(:), h0, w0, l0, channel, timepoint, dh, dw, dl);
    case 'eTypeUInt16',
        dataset.SetDataSubVolumeAs1DArrayShorts(stack(:), h0, w0, l0, channel, timepoint, dh, dw, dl);
    case 'eTypeFloat',
        dataset.SetDataSubVolumeAs1DArrayFloats(stack(:), h0, w0, l0, channel, timepoint, dh, dw, dl);
    otherwise,
        error('imarissetsubvolume: Bad value for dataset.GetType().');
end

end

