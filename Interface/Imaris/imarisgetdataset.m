function dataset = imarisgetdataset(varargin)
%
% dataset = imarisgetdataset()
%
% description: 
%   returns the current active Imaris IDataSet
% 
% output:
%   dataset     current Imaris IDataSet
%
% See also: imarissetdatasset


imaris = imarisvarargin(varargin);
dataset = imaris.GetDataSet;

end

