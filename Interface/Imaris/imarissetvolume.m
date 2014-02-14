function dataset = imarissetvolume(varargin)
%
% dataset = imarissetvolume(stack, channel, timepoint)
% dataset = imarissetvolume(object, ...)
% dataset = imarissetvolume(imaris, ...)
%
% description:
%    tries to replace volume in the active IDataSet of Imaris
%
% input:
%    stack          image stack for volumetric data
%    channel        (optional) color channel (0)
%    timepoint      (optional) timepoint to put data (0)
%
%    object         use this IDataSet object
%    imaris         Imaris application instance
%
% output:
%    dataset        IDataSet reference
%
% See also: imarisset

[imaris, varargin, nargin] = imarisvarargin(varargin);

if nargin < 1 || nargin > 4
   error('imarissetvolume: expect 1-4 input arguments!');
end

if isimaristype(imaris, varargin{1}, 'DataSet')
   dataset = varargin{1};
   varargin = varargin(2:end);
   nargin = nargin - 1;
   if nargin < 1
      error('imarissetvolume: expect at least volumetirc data as input paramter!');
   end
elseif isempty(varargin{1})
   dataset = [];
else
   dataset = imaris.GetDataSet();
end

stack = varargin{1};

% create new data set if needed or required via []
if isempty(dataset)
    
    sz = size(stack);
    if numel(sz) == 2
        sz = [sz 1];
    end
    dataset = imarissetdataset(imaris, class(stack), sz(1), sz(2), sz(3), 1, 1);

end

% transform stack if not using hwl coordinates
% stack = permute(stack, [2 1 3]); 


if nargin < 2
   channel = 0;
else
   channel = varargin{2};
end
%channel = channel - 1;
if channel > (dataset.GetSizeC() - 1)
   error('imarissetsubvolume: color index is out of bounds.');
end

if nargin < 3
   timepoint = 0;
else
   timepoint = varargin{3};
end
%timepoint = timepoint - 1;
if timepoint > (dataset.GetSizeT() - 1)
   error('imarissetsubvolume: time index is out of bounds.');
end


% Get the dataset class
switch char(dataset.GetType())
    case 'eTypeUInt8',  datatype = 'uint8';
    case 'eTypeUInt16', datatype = 'uint16';
    case 'eTypeFloat',  datatype = 'single';
    otherwise,
        error('imarissetvolume: bad value for dataset.GetType().');
end

% Check that the input and output datatypes match
if ~isa(stack, datatype)
    warning('imarissetvolume: data type mismatch. %s != %s', class(stack), datatype);
    stack = cast(stack, datatype);
end
    

% Check that the size matches
outsizes = imarisgetsize(imaris, dataset);
if ismatrix(stack)
    sizes = [size(stack) 1];
else
    sizes = size(stack);
end
if any(sizes(1 : 3) ~= outsizes(1 : 3))
    disp(sizes);
    disp(outsizes);
    error('imarissetvolume: size mismatch.');
end

% Set the stack
switch char(dataset.GetType())
    case 'eTypeUInt8',   
        dataset.SetDataVolumeAs1DArrayBytes(stack(:), channel, timepoint);
    case 'eTypeUInt16',
        dataset.SetDataVolumeAs1DArrayShorts(stack(:), channel, timepoint);
    case 'eTypeFloat',
        dataset.SetDataVolumeAs1DArrayFloats(stack(:), channel, timepoint);
    otherwise,
        error('imarissetvolume: Bad value for dataset.GetType().');
end

end

