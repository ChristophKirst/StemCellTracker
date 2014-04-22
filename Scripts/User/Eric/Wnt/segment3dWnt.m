function [imgseg, imgstats] = segment3dWnt(img, verbose)
% 3d segmentation for Wnt data

%% Init

if false
   %%
   close all
   clear all
   clc
   
   set(0, 'DefaultFigurePosition', [1   705   560   420]);
   
   initialize()
   
   ijinitialize();
end

if nargin < 2
   verbose = false;
end

%% get some test data

if false
   %%
   filename = '/home/ckirst/Science/Simulation/Matlab/StemCell3D/Test/Images/Develop/Wnt/wnt_clone8_again_feb09.lif';
   %lifdata = imread_bf();
   
   lifdata = imread_bf(filename, struct('series', 1, 'time', [3, 4], 'channel', [1, 2]));
   
   img = lifdata(:,:,:,1,1);
else
   filename = '';
end

%% initialize / prefiltering

imgd = double(img);
imgd = mat2gray(imgd);

%imglogvals = log2(imgd(:)+eps);
%imglogvals(imglogvals < -5) = -5;
%imglogvals(imglogvals > 0) = 0;

if verbose
   figure(1)
   clf
   set(gcf, 'Name', ['Raw Stack: ' filename ' channel: 1']);
   implot3d(imgd);
   
   %figure(2)
   %subplot(1,2,1);
   %hist(imglogvals, 256)
   %subplot(1,2,2);
   %hist(imgd(:), 256);
end


%% filter 
param.filter.median.ksize = [5, 5, 3]; %3;
imgmed = mat2gray(medianFilter(imgd, param.filter.median.ksize));
%imgmed = meanShiftFilter(imgmed, 4, 0.1);


imgmedf = imgmed;
imgmedf(imgmed > 0.5) = 0.5;
imgmedf = mat2gray(imgmedf);

if verbose
   %%
   figure(3)
   clf
   set(gcf, 'Name', ['Filtered Stack: ' filename ' channel: 1']);
   implot3d(imgmed);
   
   figure(4)
   clf
   set(gcf, 'Name', ['Filtered Stack: ' filename ' channel: 1']);
   implottiling(imgmedf);
   
   
   %figure(5)
   %hist(imgmed(:), 256)
end





%% thresholding / masking 

%th = 2^thresholdMixtureOfGaussians(imglogvals, 0.9);
th = thresholdMixtureOfGaussians(imgmedf, 0.5);
%th = thresholdEntropy(imgmed);
%th = thresholdMutualEntropy(imgmed);
%th = thresholdOtsu(imgmed);


imgth = imgmedf;
imgth(imgth < th) = 0;

imgmask = imgth > 0;
% remove small fragments and small isolated minima with imopen
imgmask = imopen(imgmask, strel('disk', 3));


if verbose
   %%
   figure(10)
   clf
   set(gcf, 'Name', ['Thresholded Stack: ' filename ' channel: 1']);
   implot3d(imgth);
   
   %figure(5)
   %clf
   %set(gcf, 'Name', ['Mask: ' filename ' channel: 1']);
   %imsurfaceplot3d(imgmask);
end


%% Filtering / Seeding

% gaussian smoothing
%imgdg = gaussianFilter(imgd,3,10);

% median filter  ?? also computed above
% imgf = medianFilter(imgdg, 3);

% mean shift 
%imgf = meanShiftFilter(imgd, 3, 0.1);

% disk [ksize or outer-box, width of ring, wt-inner disk, wt-ring]
%imgf = diskFilter(imgmed, [20, 20, 4], 1, 1, -1);

% Laplacian of Gaussians

zsize = size(imgmedf,3);
imgf = mat2gray(imgmedf);
for z = zsize:-1:1
   imgf(:,:,z) = imgf(:,:,z) - 0.0 * mat2gray(imgradient(imgf(:,:,z)));
end

%param.filter.logsize = [20, 20, 3];
%imgf = logFilter(max(imgf(:)) - imgf, param.filter.logsize);

imgf = sphereFilter(imgf, [20,20,4]);

%%%
% find maxima using h-max detection 
param.filter.hmax = 0.01;  %0.02;
%imgmax = imregionalmax(imgf);
imgmax = imextendedmax(mat2gray(imgf), param.filter.hmax);

if verbose
   %%
   figure(21)
   clf
   implottiling(imoverlaylabel(mat2gray(imgmedf), bwlabeln(imgmax)));
   
   figure(22)
   clf
   implottiling(imoverlaylabel(mat2gray(imgf), bwlabeln(imgmax)));
   
end


%% Watershed Segmentation on Image

%dilating the maxima can help
imgmaxd = imdilate(imgmax, strel('disk', 1));
imgmaxd = imfill(imgmaxd, 'holes');
%imgmaxd = imgmax;

% watershed
imgmin = imimposemin(max(imgmedf(:)) - imgmedf, imgmaxd);
imgws = watershed(imgmin);
imgseg = immask(imgws, imgmask);

%propagation
%imgseg = segmentByPropagation(imgmedf, bwlabeln(imgmaxd), imgmask);



%% Clean up segmentation and alternative diagnositcs
% 

param = setParameter('volume.min',    50,...     % minimal volume to keep (0)
                     'volume.max',    inf,...    % maximal volume to keep (Inf)
                     'intensity.min', -inf, ...  % minimal mean intensity to keep (-Inf)
                     'intensity.max', inf, ...   % maximal mean intensity to keep (Inf)
                     'boundaries',    false, ... % clear objects on x,y boundaries (false)
                     'fillholes',     true,...   % fill holes in each z slice after processing segments (true)
                     'relabel',       true);    % relabel from 1:nlabelnew (true)

[imgseg, stats] = postProcessSegments(imgseg, param);

if verbose
   %%
   figure(11)
   clf
   implot3d(imgpost);
   %imsurfaceplot3d(imgseg);
end

%%
if true
   ijplot3d(imcolorize(imgseg) + gray2rgb(mat2gray(img)), 'PixelDepth', 5);
end


if nargout > 1
   imgstats = stats;
end


end

