function yest = smoothLowessPeriodic(data, width, iterations, period)
%
%  out = smoothLowessPeriodic(data, width, iterations)
%
% description:
%    applies robust LOWESS smoothing filter to data assuming a periodic signal
%
% input: 
%    data        the 1d data to smooth
%    width       width for smoothing
%    iterations  number of iterations
%    period      time period (optional), xdata(end) is used otherwise
%
% output:
%    yest        the smoothed data
%
% reference:
%     William S. Cleveland, Robust locally weighted regression and smoothing scatterplots
%     Journal of the American Statistical Association, 1979, 
%     William S. Cleveland and Susan J. Devlin, Locally weighted regression: An approach 
%     to regression analysis by local fitting, Journal of the American statistical Association, 1988

if size(data, 1) == 1
   xdata = 1:numel(data);
   ydata = data;
elseif ndims(data) == 2 && size(data,2) ==1 %#ok<ISMAT>
   xdata = 1:numel(data);
   ydata = data';
elseif size(data,2) == 2
   xdata = data(:,1)';
   ydata = data(:,2)';
else
   xdata = data(1,:);
   ydata = data(2,:);
end

if nargin < 2
   width = ceil(2/3 * length(xdata));
end
if width >= length(xdata)
   error('smoothLowessPeriodic: width of smoothing bigger that nubmer of data points!')
end
if nargin < 3
   iterations = 1;
end
if nargin < 4
   period = xdata(end);
end

xdata = [xdata((end - width) : end) - period, xdata, xdata(1:width+1) + period] + period;
ydata = [ydata((end - width) : end),          ydata, ydata(1:width+1)];

yest = smoothLowess([xdata; ydata], width, iterations);

yest = yest(width+2:end-width-1);

end
