function [pks, img] = detectPeaksByClosing(img, varargin)
%
% pks = detectPeaksByClosing(img, varargin)
%
% description:
%    find peaks in image using the specified parameter
%
% input:
%   img     image
%   param   (optional) parameter struct with entries
%           .filtersize Gaussian filter size ([] = no gaussian filtering)
%           .strel      structural element for closing
%           .scale      rescaling factor
%           .threshold  threshold after closing ([] = thresholdFirstMin)
%
% output:
%   pks      array of pixel coordinates (coordinates of each peak as one row vector)

param = parseParameter(varargin);

%
rs = getParameter(param, 'scale', []);
if ~isempty(rs)
   img = imresize(img, rs);
end

% filter 
fs = getParameter(param, 'filtersize', []);
if ~isempty(fs)
   img = filterGaussian(img, fs);
end

% morphological closing
se = getParameter(param, 'strel', strel('disk', 20));
if isnumeric(se)
   se = strel('disk', se);
end
img = imclose(img, se);

{min(img(:)), max(img(:))}

th = getParameter(param, 'threshold', []);
if isempty(th)
   th = thresholdFirstMin(mat2gray(img));
end

img = img > th;

[idx, idy] = find(img);

pks = [idx'; idy'];

% % plot 
% figure
% implot(img)
% hold on
% plot(idx, idy, '*', 'Color', 'r')


if ~isempty(rs)
   pks = (pks-1) * 1/rs  + 1;
end

fprintf('detectPeaksByClosing: found %g peaks\n', size(pks, 2));

end