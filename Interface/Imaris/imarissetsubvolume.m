function dataset = imarissetsubvolume(varargin)
%
% stack = imarissetsubvolume(stack, p0, q0, l0, channel, timepoint)
% stack = imarissetsubvolume(stack, [p0, q0, l0], channel, timepoint)
% stack = imarissetsubvolume(object, ...)
% stack = imarissetsubvolume(imaris, ...)
%
% description:
%    tries to replace a sub volume in the active IDataSet of Imaris
%
% input:
%    pql0           position in pixel coordinates
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


pql0 = varargin{2};
if length(pql0) == 1
   if nargin < 4 
      error('imarissetsubvolume: expect 3 parameters: p0, q0, l0');
   end
   pql0 = [varargin{2:4}];
   varargin = varargin(5:end);
   nargin = length(varargin);
   
elseif length(pql0) == 3
   
   varargin = varargin(3:end);
   nargin = length(varargin);
   
else 
   error('imarissetsubvolume: input syntax error.');
   
end

% matlab first index = 1, Imaris first index = 0
pql0 = pql0 -1;
p0 = pql0(1); q0 = pql0(2); l0 = pql0(3);

if p0 < 0 || p0 > dataset.GetSizeX() - 1
   error('imarisget: starting position p0 out of bounds.');
end
p0 = uint32(p0);

if q0 < 0 || q0 > dataset.GetSizeY() - 1
   error('imarissetsubvolume: starting position q0 out of bounds.');
end
q0 = uint32(q0);

if l0 < 0 || l0 > dataset.GetSizeZ() - 1
   error('imarissetsubvolume: starting position l0 out of bounds.');
end
l0 = uint32(l0);

 
[dp, dq, dl] = size(stack);

dsize = imarisgetsize(imaris, dataset);

if any(pql0 + [dp,dq,dl] > dsize -1)
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
        dataset.SetDataSubVolumeAs1DArrayBytes(stack(:), p0, q0, l0, channel, timepoint, dp, dq, dl);
    case 'eTypeUInt16',
        dataset.SetDataSubVolumeAs1DArrayShorts(stack(:), p0, q0, l0, channel, timepoint, dp, dq, dl);
    case 'eTypeFloat',
        dataset.SetDataSubVolumeAs1DArrayFloats(stack(:), p0, q0, l0, channel, timepoint, dp, dq, dl);
    otherwise,
        error('imarissetsubvolume: Bad value for dataset.GetType().');
end

end

