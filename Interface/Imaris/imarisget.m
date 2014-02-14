function data = imarisget(varargin)
%
% data = imarisget('Volume', channel, timepoint, iDataSet)
% data = imarisget(imarisApplication, 'Volume', channel, timepoint, iDataSet)
%
% description:
%    sets a stack in Imaris in color channel and timepoint to stack
%    or setst a surface object defined by coordinates X and facets K
%
% input:
%   channel  : channel number 
%   timepoint: timepoint number
%   iDataSet : (optional) get the data volume from the passed IDataset
%               object instead of current one; if omitted, current dataset
%               (i.e. .mImarisApplication.GetDataSet()) will be used.
%               This is useful for instance when masking channels.
%
% output:
%   data        stack
%
% TODO: surfaces, spots, other usefull objects
%
% See also: imarisput

[mImarisApplication, var] = imarisvarargin(varargin);
nargin = length(varargin);
varargin = var{:};


if nargin < 1 
   type = 'Volume';
else
   type = varargin{1};
   varargin = varargin(2:end);
   nargin = nargin -1;
end


typeValues = {'Volume', 'SubVolume'};

if ~ismember(type, typeValues)
   error('imarisget: data type is not supported!')
end


switch type
   case 'Volume'
   
      % Initialize stack
      data = [];

      if nargin < 1
         channel = 0;
      else
         channel = varargin{1};
      end
      if nargin < 2
         timepoint = 0;
      else
         timepoint = varargin{2};
      end
      if nargin < 3     
          iDataSet = mImarisApplication.GetDataSet();
      else
          iDataSet = varargin{3};
          % Is the passed dataset a valid DataSet?
          if ~mImarisApplication.GetFactory().IsDataSet(iDataSet)
              error('Invalid IDataSet object.');
          end
      end

      % Check whether we have some voxels in the dataset
      if isempty(iDataSet) || iDataSet.GetSizeX() == 0
          return
      end

      % Convert channel and timepoint to 0-based indexing
      %channel = channel - this.mIndexingStart;
      %timepoint = timepoint - this.mIndexingStart;

      % Check that the requested channel and timepoint exist
      if channel > (iDataSet.GetSizeC() - 1)
          error('The requested channel index is out of bounds.');
      end
      if timepoint > (iDataSet.GetSizeT() - 1)
          error('The requested time index is out of bounds.');
      end

      % Get the dataset class
      switch char(iDataSet.GetType())
          case 'eTypeUInt8',   datatype = 'uint8';
          case 'eTypeUInt16',  datatype = 'uint16';
          case 'eTypeFloat',   datatype = 'single';
          otherwise,
              error('Bad value for IDataSet::GetType().');
      end

      % Allocate memory
      data = zeros([iDataSet.GetSizeX(), iDataSet.GetSizeY(), ...
          iDataSet.GetSizeZ()], datatype);

      % Get the stack
      switch char(iDataSet.GetType())
          case 'eTypeUInt8',   
              % Java does not have unsigned ints
              %channel
              %timepoint
              arr = iDataSet.GetDataVolumeAs1DArrayBytes(channel, timepoint);
              data(:) = typecast(arr, 'uint8');
          case 'eTypeUInt16',
              % Java does not have unsigned ints
              arr = iDataSet.GetDataVolumeAs1DArrayShorts(channel, timepoint);
              data(:) = typecast(arr, 'uint16');
          case 'eTypeFloat',
              data(:) = ...
                  iDataSet.GetDataVolumeAs1DArrayFloats(channel, timepoint);
          otherwise,
              error('Bad value for iDataSet.GetType().');
      end
      
      %data = permute(data,[2,1,3]);
      
    
    case 'SubVolume'   % (x0, y0, z0, dx, dy, dz, channel, timepoint, iDataSet)
       
       
       data = [];
       
       if nargin < 7 || nargin > 9
          error('imarisget: expects 7 to 9 input arguments!');
       end
       
       
       if nargin < 9
          iDataSet = mImarisApplication.GetDataSet();
       else
          iDataSet = varargin{9};
          if ~mImarisApplication.GetFactory().IsDataSet(iDataSet)
             error('imarisget: invalid IDataset object.');
          end
       end
       
       % Check whether we have some voxels in the dataset
       if isempty(iDataSet) || iDataSet.GetSizeX() == 0
          return
       end
       
       % Convert all dimensions to 0-based indexing
       indexOffset = 0;
       x0 = varargin{1} - indexOffset;
       if x0 < 0 || x0 > iDataSet.GetSizeX() - 1
          error('imarisget: starting position x0 out of bounds.');
       end
       x0 = uint32(x0);
       dx = varargin{4};

       y0 = varargin{2} - indexOffset;
       if y0 < 0 || y0 > iDataSet.GetSizeY() - 1
          error('imarisget: starting position y0 out of bounds.');
       end
       y0 = uint32(y0);
       dy = varargin{5};
       
       z0 = varargin{3} - indexOffset;
       if z0 < 0 || z0 > iDataSet.GetSizeZ() - 1
          error('imarisget: starting position z0 out of bounds.');
       end
       z0 = uint32(z0);
       dz = varargin{6};
       
       channel = varargin{7} - indexOffset;
       if channel < 0 || channel > iDataSet.GetSizeC() - 1
          error('imarisget: channel index out of bounds.');
       end
       channel = uint32(channel);
       
       if nargin < 8
          timepoint = 0;
       else
          timepoint = varargin{8};
       end
       if timepoint < 0 || timepoint > iDataSet.GetSizeT() - 1
          error('imarisget: timepoint index out of bounds.');
       end
       timepoint = uint32(timepoint);
       
       % Check that we are within bounds
       if x0 + dx > iDataSet.GetSizeX()
          error('imarisget: x range is out of bounds.');
       end
       
       if y0 + dy > iDataSet.GetSizeY()
          error('imarisget: y range is out of bounds.');
       end
       
       if z0 + dz > iDataSet.GetSizeZ()
          error('imarisget: z range is out of bounds.');
       end
       
       % Get the dataset class
       switch char(iDataSet.GetType())
          case 'eTypeUInt8',   datatype = 'uint8';
          case 'eTypeUInt16',  datatype = 'uint16';
          case 'eTypeFloat',   datatype = 'single';
          otherwise,
             error('imarisget: Bad value for IDataSet::GetType().');
       end
       
       % Allocate memory
       data = zeros([dx, dy, dz], datatype);
       
       % Get the stack
       switch char(iDataSet.GetType())
          case 'eTypeUInt8',
             % Java does not have unsigned ints
             arr = iDataSet.GetDataSubVolumeAs1DArrayBytes(x0, y0, z0, ...
                channel, timepoint, dx, dy, dz);
             data(:) = typecast(arr, 'uint8');
          case 'eTypeUInt16',
             % Java does not have unsigned ints
             arr = iDataSet.GetDataSubVolumeAs1DArrayShorts(x0, y0, z0, ...
                channel, timepoint, dx, dy, dz);
             data(:) = typecast(arr, 'uint16');
          case 'eTypeFloat',
             data(:) = ...
                iDataSet.GetDataSubVolumeAs1DArrayFloats(x0, y0, z0, ...
                channel, timepoint, dx, dy, dz);
          otherwise,
             error('imarisget: Bad value for iDataSet.GetType().');
       end
       
       %data = permute(data,[2,1,3]);
end

      
      
%       
%     case 'SubVolume'
%       % Initialize stack
%       data = [];
%       
%       if nargin < 1
%           x0 = 0;
%       else
%           x0 = varargin{1};
%       end
%       
%       
%       
%       
%       if nargin < 1
%          channel = 0;
%       else
%          channel = varargin{1};
%       end
%       if nargin < 2
%          timepoint = 0;
%       else
%          timepoint = varargin{2};
%       end
%       if nargin < 3     
%           iDataSet = mImarisApplication.GetDataSet();
%       else
%           iDataSet = varargin{3};
%           % Is the passed dataset a valid DataSet?
%           if ~mImarisApplication.GetFactory().IsDataSet(iDataSet)
%               error('Invalid IDataSet object.');
%           end
%       end
% 
%       % Check whether we have some voxels in the dataset
%       if isempty(iDataSet) || iDataSet.GetSizeX() == 0
%           return
%       end
% 
%       % Convert channel and timepoint to 0-based indexing
%       %channel = channel - this.mIndexingStart;
%       %timepoint = timepoint - this.mIndexingStart;
% 
%       % Check that the requested channel and timepoint exist
%       if channel > (iDataSet.GetSizeC() - 1)
%           error('The requested channel index is out of bounds.');
%       end
%       if timepoint > (iDataSet.GetSizeT() - 1)
%           error('The requested time index is out of bounds.');
%       end
% 
%       % Get the dataset class
%       switch char(iDataSet.GetType())
%           case 'eTypeUInt8',   datatype = 'uint8';
%           case 'eTypeUInt16',  datatype = 'uint16';
%           case 'eTypeFloat',   datatype = 'single';
%           otherwise,
%               error('Bad value for IDataSet::GetType().');
%       end
% 
%       % Allocate memory
%       data = zeros([iDataSet.GetSizeX(), iDataSet.GetSizeY(), ...
%           iDataSet.GetSizeZ()], datatype);
% 
%       % Get the stack
%       switch char(iDataSet.GetType())
%           case 'eTypeUInt8',   
%               % Java does not have unsigned ints
%               arr = iDataSet.GetDataVolumeAs1DArrayBytes(channel, timepoint);
%               data(:) = typecast(arr, 'uint8');
%           case 'eTypeUInt16',
%               % Java does not have unsigned ints
%               arr = iDataSet.GetDataVolumeAs1DArrayShorts(channel, timepoint);
%               data(:) = typecast(arr, 'uint16');
%           case 'eTypeFloat',
%               data(:) = ...
%                   iDataSet.GetDataVolumeAs1DArrayFloats(channel, timepoint);
%           otherwise,
%               error('Bad value for iDataSet.GetType().');
%       end
%       
%       data = permute(data,[2,1,3]);
      
end