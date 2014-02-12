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

[mImarisApplication, varargin] = imarisvarargin(varargin);
nargin = length(varargin);

if nargin < 1 
   type = 'Volume';
else
   type = varargin{1};
   varargin = {varargin(2:end)};
end


typeValues = {'Volume'};

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
      
      data = permute(data,[2,1,3]);
      
end