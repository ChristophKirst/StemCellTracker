function threshold = thresholdFirstMin(image, varargin)
%
% threshold = thresholdFirstMin(image, param)
%
% description:
%    finds a threshold by searching first minima of the histogram
%
% input:
%    image       the image to calculate threshold for
%    param       paramter struct with entry
%                .nbins     number of bins in the histogram
%                .delta     height of minimum
%
% output:
%    threshold   calculated threshold
%

param =parseParameter(varargin);

image = double(image);

%histrogram

[N, X] = hist(image(:), getParameter(param, 'nbins', 50));

tid = findTroughs(N, getParameter(param, 'delta', 1/100 * N), 1);

if isempty(tid)
   threshold = X(end);
else
   threshold = X(tid);
end
   
end











