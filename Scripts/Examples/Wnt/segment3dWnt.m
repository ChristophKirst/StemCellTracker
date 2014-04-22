function imgpost = segment3dWnt(img, verbose)
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
   
   bfinitialize;
   
   verbose = true;
end

if nargin < 2
   verbose = false;
end

%% get some test data

if false
   %%
   filename = '/home/ckirst/Science/Projects/StemCells/Experiment/Other/Wnt/wnt_clone8_again_feb09.lif';
   %lifdata = imread_bf();
   
   lifdata = imread_bf(filename, struct('series', 1, 'time', [3, 4], 'channel', [1, 2])); 
   img = lifdata(:,:,:,1,1);
   
   
   %add a zero image on top
   %img(:,:,end+1) = zeros(size(img(:,:,end)));
   
else
   filename = '';
end

%% initialize / prefiltering

img = double(img);
img = mat2gray(img);

%imglogvals = log2(imgd(:)+eps);
%imglogvals(imglogvals < -5) = -5;
%imglogvals(imglogvals > 0) = 0;

if verbose
   figure(1)
   clf
   set(gcf, 'Name', ['Raw Stack: ' filename ' channel: 1']);
   implot3d(imresample(img, [3, 3, 1]));
   %ijplot3d(imgd, 'PixelDepth', 8)
   
   %figure(2)
   %subplot(1,2,1);
   %hist(imglogvals, 256)
   %subplot(1,2,2);
   %hist(imgd(:), 256);
end


%% filter 
imgf = mat2gray(medianFilter(img, [5, 5, 3], 'replicate'));
%imgmed = meanShiftFilter(imgmed, 4, 0.1);

imgf(imgf > 0.5) = 0.5;
imgf = mat2gray(imgf);

if verbose   
   figure(3)
   clf
   set(gcf, 'Name', ['Filtered Stack: ' filename ' channel: 1']);
   implottiling(imgf);
  
   %figure(4)
   %hist(imgmed(:), 256)
end


%% thresholding / masking 

%th = 2^thresholdMixtureOfGaussians(imglogvals, 0.9);
th = thresholdMixtureOfGaussians(imgf, 0.5);
%th = thresholdEntropy(imgmed);
%th = thresholdMutualEntropy(imgmed);
%th = thresholdOtsu(imgmed);


imgth = imgf;
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
   
   figure(5)
   clf
   set(gcf, 'Name', ['Mask: ' filename ' channel: 1']);
   implot3d(imgmask);
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

zsize = size(imgf,3);
imgs = mat2gray(imgf);
for z = zsize:-1:1
   imgs(:,:,z) = imgs(:,:,z) - 0.1 * mat2gray(imgradient(imgth(:,:,z)));
end

%imgs = logFilter(max(imgs(:)) - imgs, [17,17,3], [], 'replicate');
%imgs = sphereFilter(imgs, [22,22,4], [], 'replicate');
imgs = diskFilter(imgs, [10,10,3]);

%%%
% find maxima using h-max detection 
%imgmax = imregionalmax(imgs);
imgmax = imextendedmax(mat2gray(imgs), 0.05);

%dilating the maxima can help
%imgmax = imdilate(imgmax, strel('disk', 1));
%imgmax = imfill(imgmax, 'holes');
%imgmaxd = imgmax;

if verbose
   %%
   figure(21); clf
   set(gcf, 'Name', ['Seeding: ' filename ' channel: 1']);
   implottiling(imoverlay(mat2gray(imgs), bwlabeln(imgmax)));
   
   figure(22); clf
   set(gcf, 'Name', ['Seeding / Thresholded Image: ' filename ' channel: 1']);
   implottiling(imoverlaylabel(mat2gray(imgth), bwlabeln(imgmax)));
   
end


%% Segmentation on Image

% watershed
imgmin = imimposemin(max(imgs(:)) - imgs, imgmax);
imgws = watershed(imgmin);

if verbose 
   %%
   figure(42)
   implot3d(imcolorize(imgws))

   %%
   figure(43)
   imgseg = immask(imgws, imgmask);
   implot3d(imcolorize(imgseg))

   %%
   figure(44)
   implotlabeloutline(imgf, imgseg)
end

%% Postprocessing

param = setParameter('volume.min',    50,...     % minimal volume to keep (0)
                     'volume.max',    inf,...    % maximal volume to keep (Inf)
                     'intensity.min', -inf, ...  % minimal mean intensity to keep (-Inf)
                     'intensity.max', inf, ...   % maximal mean intensity to keep (Inf)
                     'boundaries',    false, ... % clear objects on x,y boundaries (false)
                     'fillholes',     true,...   % fill holes in each z slice after processing segments (true)
                     'relabel',       true);    % relabel from 1:nlabelnew (true)

[imgpost, stats] = postProcessSegments(imgseg, param);

if verbose
   
   figure(47); clf
   implottiling(imcolorize(imgpost))
end


%%
if verbose
   ijplot3d(imcolorize(imgseg) + gray2rgb(mat2gray(img)), 'PixelDepth', 4);
end


end

