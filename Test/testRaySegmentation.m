%% Test Segmentation by Rays

clear all
close all
clc


%% load some image data

img = imread('./Test/Images/raw.tif');
[imggrad, imggraddir] = imgradient(gaussianFilter(img,3));
imggrad = mat2gray(imggrad);

% plot
% figure(1);
% clf
% ax(1) = imsubplot(1,2,1);
% imshow(img)
% title('raw')
% ax(2) = imsubplot(1,2,2);
% imshow(imggrad);
% title('gradient')
% linkaxes(ax, 'xy');


% extract smaller area for testing and filter noise

crop = [166.5 207.5 194.5 191.5];
imgd = imcrop(img, crop);
imgdgrad =  imcrop(imggrad, crop);

imgd = mat2gray(medianFilter(imgd,3));

%imgd = gaussianFilter(imgd,3, 100);
%imgdgrad = gaussianFilter(imgdgrad,3, 100);

% plot
figure(1)
clf
ax(1) = imsubplot(1,2,1);
imshow(imgd)
title('raw')
ax(2) = imsubplot(1,2,2);
imshow(imgdgrad);
title('gradient')
linkaxes(ax, 'xy');

%% detect mask

imglogvals = log2(imgd(:)+eps);
imglogvals(imglogvals < -15) = -15;
imglogvals(imglogvals > 0) = 0;
%th = 2^thresholdMixtureOfGaussians(imglogvals, 0.9);
th = 2^-5;

imgth = imgd;
imgth(imgth < th) = 0;

imgmask = imgth > 0;
imgmask = imopen(imgmask, strel('disk', 1));
imgmask = imclose(imgmask, strel('disk', 1));


figure(2)
clf
imsubplot(1,3,1)
imshow(imgth);
imsubplot(1,3,2)
imshow(imoverlay(imgd, imgmask))
subplot(1,3,3)
hist(imglogvals, 256)



%% get some points

npoints = 0;

clear pt
figure(1)
for p=npoints:-1:1
   pt(p) = impoint();         % retunrs spatial coordinates !!!
   r0(:,p) = round(pt(p).getPosition())';
   r0(:,p) = r0(2:-1:1, p);   % to indices
end

%r0 = [ 117    31    97;  48   165    75];

r0 = [79    66   120   158   152    99    85    34    79    31;
      12   132   150   123    52    78    51    57    39   164];

r0 = r0(:,:);   
   
seeds = index2mask(size(imgd), r0(:,:));

% plot
figure(1)
clf
ax(1) = imsubplot(1,2,1);
imshow(imoverlay(imgd, seeds));
title('raw')
ax(2) = imsubplot(1,2,2);
imshow(imoverlay(imgdgrad, seeds));
title('gradient')
linkaxes(ax, 'xy');


%% test findRayRadius.m on seeds


seeds = index2mask(size(imgd), r0);

param.nrays = 25;
param.cutoff.radius = 20;
param.center_averaging = 2;


param.threshold.background = 0;
%param.threshold.background = 0;

param.threshold.relative_change = 1.0;        % stop propagation if reltative change is larger (0.5, Inf = none)
%param.threshold.relative_change = Inf;

%param.threshold.absolute_change = 1000000;       % stop propagation if absolute change is larger (Inf = none)
param.threshold.absolute_change = Inf;

param.threshold.relative_gradient_peak = Inf; % minimal peak height in gradient profile (Inf = none)
%param.threshold.relative_gradient_peak = Inf

param.threshold.absolute_gradient_peak = 0.15; % minimal peak height in gradient profile (Inf = none)
%param.threshold.absolute_gradient_peak = Inf;

param.trough.drop = 0.25;                     % minimal trough falloff in itensity profile (Inf = do not use trough method)
%param.trough.drop = Inf; 
param.trough.max_rise = 0.5;                  % stop trough finding if profile rises more than this
param.trough.rise = 0.5;                      % minimal subsequent trough increase in itensity profile
param.trough.noise = 0.2;                     % if no increase reset trough minima by ignoring fluctuations less than this         



param.smooth = 6;                             % width of LOWESS filter for smoothing


param.plot = true;

findRayRadius(imgth, imgdgrad, seeds, param);




%% detect a bunch of maxima

% gaussian smoothing
%imgdg = gaussianFilter(imgd,3,10);
imgdg = imgd;


% median filter
%imgf = medianFilter(imgdg, 3);

% mean shift 
%imgf = meanShiftFilter(imgd, 3, 0.1);

% disk
%imgf = diskFilter(imgmed, 3, 4, 1, -1);

% Laplcaina of Gaussians
param.filter.logsize = [15, 15];
imgf = logFilter(max(imgdg(:)) - imgdg, param.filter.logsize);


% find maxima using h-max detection 
param.filter.hmax = 0.005;
%imgmax = imregionalmax(imgf);
imgmax = imextendedmax(mat2gray(imgf), param.filter.hmax);
imgmax = imgmax .* imgmask;

%%% Combine nearby points and shrink to single points
imgmax = imdilate(imgmax, strel('disk', 3));
imgmax = bwmorph(imgmax,'shrink',inf);


% plot the results
figure(10)
imshow([imoverlay(imgd, imgmax) imoverlay(mat2gray(imgf), imgmax)])



%% test findRayRadius.m on found seeds 

seeds = imgmax;

param.nrays = 25;
param.cutoff.radius = 20;
param.center_averaging = 2;


param.threshold.background = 0;
%param.threshold.background = 0;

param.threshold.relative_change = 1.0;        % stop propagation if reltative change is larger (0.5, Inf = none)
%param.threshold.relative_change = Inf;

%param.threshold.absolute_change = 1000000;       % stop propagation if absolute change is larger (Inf = none)
param.threshold.absolute_change = Inf;

param.threshold.relative_gradient_peak = Inf; % minimal peak height in gradient profile (Inf = none)
%param.threshold.relative_gradient_peak = Inf

param.threshold.absolute_gradient_peak = 0.15; % minimal peak height in gradient profile (Inf = none)
%param.threshold.absolute_gradient_peak = Inf;

param.trough.drop = 0.25;                     % minimal trough falloff in itensity profile (Inf = do not use trough method)
%param.trough.drop = Inf; 
param.trough.max_rise = 0.5;                  % stop trough finding if profile rises more than this
param.trough.rise = 0.5;                      % minimal subsequent trough increase in itensity profile
param.trough.noise = 0.2;                     % if no increase reset trough minima by ignoring fluctuations less than this         



param.smooth = 6;                             % width of LOWESS filter for smoothing


param.plot = true;

findRayRadius(imgth, imgdgrad, seeds, param);

% oversegmentation is corrected here automatically,
% the algorithm runs slower than others though


 

%% todo: statistics 

mask = poly.createMask;
%figure(2)
%clf
%imshow(mask)

stats.area.pos = find(mask);
stats.area.values = imgz(stats.area.pos);
stats.area.size = sum(sum(mask));
stats.area.intensity.mean = mean(double(stats.area.values));
stats.area.intensity.std = std(double(stats.area.values), 1);


stats.border.pos = find(bwperim(mask));
stats.border.values = imgz(stats.border.pos);
stats.border.length = length(stats.border.pos);
stats.border.grad.values   = imgzgrad(stats.border.pos);
stats.border.grad.mean = mean(stats.border.grad.values);
stats.border.grad.std = std(double(stats.border.grad.values));


stats.pos = poly.getPosition;

pos2 = zeros(nlines,2);
pos2(1:end-1, :) = stats.pos(2:end, :);
pos2(end, :) = stats.pos(1, :);

stats.distance = sqrt(sum((stats.pos - pos2).^2,2));
