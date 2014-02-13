function varargout = imarisgetsize(varargin)
%
% varargout = imarisgetsize(varargin)
%
% 

[mImarisApplication, varargin, nargin] = imarisvarargin(varargin{:});

if isempty(mImarisApplication.GetDataSet())
    error('imarisgetsize: not valid data set');
end

if nargout <= 1
   if nargin == 1 && strcmp(varargin{1}, 'all')
      varargout{1} = [ mImarisApplication.GetDataSet.GetSizeX, ...
                       mImarisApplication.GetDataSet.GetSizeY, ...
                       mImarisApplication.GetDataSet.GetSizeZ, ...
                       mImarisApplication.GetDataSet.GetSizeC, ...
                       mImarisApplication.GetDataSet.GetSizeT ];
   else
      varargout{1} = [ mImarisApplication.GetDataSet.GetSizeX, ...
                       mImarisApplication.GetDataSet.GetSizeY, ...
                       mImarisApplication.GetDataSet.GetSizeZ ];
   end
else    
    varargout{1} = mImarisApplication.GetDataSet.GetSizeX;
    varargout{2} = mImarisApplication.GetDataSet.GetSizeY;
    varargout{3} = mImarisApplication.GetDataSet.GetSizeZ;
    varargout{4} = mImarisApplication.GetDataSet.GetSizeC;
    varargout{5} = mImarisApplication.GetDataSet.GetSizeT;
end
