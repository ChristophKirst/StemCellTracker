function imarissave(varargin)
%
% imarissave(filename, opts)
% imarissave(imaris, filename, opts)
%
% description:
%    opens a file in imaris
%
% input:
%    filename   file to save to
%    imaris     (optional) Imaris Application instance
%    opts       (optional) options ('')
%      
% usage:
%    Formats:    writer="BMPSeries"
%                Available formats:  Imaris5, Imaris3, Imaris2, SeriesAdjustable, TiffSeriesRGBA, ICS, OlympusCellR, OmeXml, BMPSeries, MovieFromSlices.
% 
% See also: imarisopen

[imaris, varargin, nargin] = imarisvarargin(varargin);

if nargin < 1
   error('imarissave: need to sepcify filename')
else
   filename = varargin{1};
end
if nargin < 2
   opts = '';
else
   opts = varargin{2};
end
   
imaris.FileSave(filename, opts);

end

