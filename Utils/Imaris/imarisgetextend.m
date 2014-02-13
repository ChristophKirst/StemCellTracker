function varargout = imarisgetextend(varargin)
%
% varargin = imarisgetextend(varargout)
%
% description:
%    get spacial extend of the image
%
% See also: imarisgetsize

mImarisApplication = imarisvarargin(varargin{:});
mDataSet = mImarisApplication.GetDataSet();

if isempty(mDataSet)
    error('imarisgetsize: not valid data set');
end

if nargout <= 1
   varargout{1} = [ 
      mDataSet.GetExtendMinX, mDataSet.GetExtendMinY, mDataSet.GetExtendMinZ; ...
      mDataSet.GetExtendMaxX,  mDataSet.GetExtendMaxY, mDataSet.GetExtendMaxZ ];
else
   varargout{1} = [mDataSet.GetExtendMinX, mDataSet.GetExtendMaxX];
   varargout{2} = [mDataSet.GetExtendMinY, mDataSet.GetExtendMaxY];
   varargout{3} = [mDataSet.GetExtendMinZ, mDataSet.GetExtendMaxZ];
end

end