function [pks, img] = detectPeaksByHmax(img, varargin)
%
% pks = detectPeaksByHmax(img, varargin)
%
% description:
%    find peaks in image usin the specified parameter
%
% input:
%   img     image
%   param   (optional) parameter struct with entries
%           .filtersize Gaussian filter size ([] = no gaussian filtering)
%           .hmax       h value for hmax transform ([] = 1/10 median(img))
%           .center     true = return center of peaks only, false = return coordinatss of all peak pixels (true) 
%
% output:
%   pks      array of pixel coordinates (coordinates of each peak as one row vector)

param = parseParameter(varargin);

% filter 
fs = getParameter(param, 'filtersize', []);
if ~isempty(fs)
   img = filterGaussian(img, fs);
end

% find extended maxima
h = getParameter(param, 'hmax', []);
if isempty(h)
   h = 1/10 * median(img(:));
end
img = imextendedmax(img, h);

% center
c = getParameter(param, 'centroid', true);
if c
   img = bwlabeln(img);
   stats = imstatistics(img, 'Centroid');
   pks = [stats.Centroid];
else
   pks = imstatistics(img, 'PixelList');
   pks = [pks.PixelList];
end

end