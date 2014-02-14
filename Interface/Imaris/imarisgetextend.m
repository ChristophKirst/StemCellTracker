function varargout = imarisgetextend(varargin)
%
% varargin = imarisgetextend(varargout)
%
% description:
%    get spacial extend of the image
%
% See also: imarisgetsize

[imaris, varargin, nargin] = imarisvarargin(varargin{:});

if nargin > 1
   if isimaristye(varargin{1}, 'DataSet')
      dataset = varargin{1};
   else
      error('imarisgetextend: passed object is not a DataSet.')
   end
else
   dataset = imaris.GetDataSet();
end

if isempty(dataset)
    error('imarisgetextend: no active data set found.');
end

if nargout == 1
   varargout{1} = [ 
      dataset.GetExtendMinX, dataset.GetExtendMinY, dataset.GetExtendMinZ; ...
      dataset.GetExtendMaxX,  dataset.GetExtendMaxY, dataset.GetExtendMaxZ ];
   
elseif nargout == 2
   varargout{1} = [dataset.GetExtendMinX, dataset.GetExtendMinY, dataset.GetExtendMinZ];
   varargout{2} = [dataset.GetExtendMaxX,  dataset.GetExtendMaxY, dataset.GetExtendMaxZ];
   
else
   varargout{1} = [dataset.GetExtendMinX, dataset.GetExtendMaxX];
   varargout{2} = [dataset.GetExtendMinY, dataset.GetExtendMaxY];
   varargout{3} = [dataset.GetExtendMinZ, dataset.GetExtendMaxZ];
end

end