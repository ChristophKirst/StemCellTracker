function varargout = imarisgetsize(varargin)
%
% isize = imarisgetsize(varargin)
% [ixsize, iysize, izsize, icsize, itsize] = imarisgetsize(varargin)
%
% description:
%    returns the size of the current IDataSet in Imaris
%
% input:
%    dataset   (optional) IDataSet
%    'all'     (optional) return full data size (x,y,z,c,t)
%
% output:
%    isize     x,y,z size or x,y,z,c,t size
%
% See also: imarisgetextend

[imaris, varargin, nargin] = imarisvarargin(varargin);

if nargin > 1
   if isimaristye(varargin{1}, 'DataSet')
      dataset = varargin{1};
   else
      error('imarisgetdatasize: passed object is not a DataSet.')
   end
else
   dataset = imaris.GetDataSet();
end

if isempty(dataset)
    error('imarisgetvoxelsize: no active data set found.');
end

if nargout == 1
   if nargin == 1 && ischar(varargin{1}) && strcmp(varargin{1}, 'all')
      varargout{1} = [ dataset.GetSizeX, ...
                       dataset.GetSizeY, ...
                       dataset.GetSizeZ, ...
                       dataset.GetSizeC, ...
                       dataset.GetSizeT ];
   else
      varargout{1} = [ dataset.GetSizeX, ...
                       dataset.GetSizeY, ...
                       dataset.GetSizeZ ];
   end
else    
    varargout{1} = dataset.GetSizeX;
    varargout{2} = dataset.GetSizeY;
    varargout{3} = dataset.GetSizeZ;
    varargout{4} = dataset.GetSizeC;
    varargout{5} = dataset.GetSizeT;
end
