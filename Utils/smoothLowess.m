function yest = smoothLowess(data, width, iterations)
%
%  out = smoothLowess(data, width_fraction, iterations)
%
% description:
%    applies robust LOWESS smoothing filter to data
%
% input: 
%    data        the 1d data to smooth
%    width       width for smoothing
%    iterations  number of iterations
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
   n = numel(data);
   xdata = 1:n;
   ydata = data;
elseif ndims(data) == 2 && size(data,2) ==1 %#ok<ISMAT>
   n = numel(data);
   xdata = 1:n;
   ydata = data';
elseif size(data,2) == 2
   n = size(data,1);
   xdata = data(:,1)';
   ydata = data(:,2)';
else
   n = size(data,2);
   xdata = data(1,:);
   ydata = data(2,:);
end

if nargin < 2
   width = ceil(2/3 * n);
end
if nargin < 3
   iterations = 1;
end


for i = n:-1:1
   so = sort(abs(xdata - xdata(i)));
   h(i) = so(width+1);
end


w = abs(bsxfun(@minus, xdata, xdata') ./ repmat(h, n, 1));
w(w > 1) = 1.0; w(w < 0) = 0.0;

w = 1 - w .* w .* w;
w = w .* w .* w;
yest = zeros(1,n);
delta = ones(1,n);
 
for iter = 1:iterations
   
   for i = 1:n 
      weights = delta .* w(i, :);
      
      weights_mul_x = weights .* xdata; 
      sum_weights = sum(weights);
      
      meanX  = sum(weights_mul_x) / sum_weights;
      meanY  = dot(weights, ydata) / sum_weights;
      meanX2 = dot(weights_mul_x, xdata) / sum_weights;
      meanXY = dot(weights_mul_x, ydata) / sum_weights;
      
      meanXmeanX = meanX * meanX;
      
      if (meanX2 == meanXmeanX)
         beta = 0.0;
      else
         beta = (meanXY - meanX * meanY) / (meanX2 - meanXmeanX);
      end
      
      alpha = meanY - beta * meanX;
      
      yest(i) = alpha + beta * xdata(i);  
   end   
   
   
   if iter == iterations
      return
   end
   
   residuals = ydata - yest;
   s = median(abs(residuals));

   delta = residuals ./ (6 * s);
   delta(delta > 1) = 1; delta(delta < -1) = -1;
   delta = 1 - delta .* delta;
   delta = delta .* delta;
end

end


%% debugging
% clc
%  
% y = [1 2 3 5 8 3 1 -5 -10 -2 0];
% ys = smoothLowess(y, 3,1);
% ysp = smoothLowessPeriodic(y, 5, 1);
% 
% figure(42)
% clf
% hold on
% plot(y, '*')
% plot(ys, '+r')
% plot(ys, 'r')
% plot(ysp, 'og')
% plot(ysp, 'g')

