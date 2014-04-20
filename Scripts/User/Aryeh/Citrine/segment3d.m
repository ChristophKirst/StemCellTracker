function [imgpost, stats] = segment3d(img, verbose, figure_offset)
%
% skript for segmenting the Citrine data
%

%% Initialize 
if false
   %%
   close all
   clear all
   clc
   set(0, 'DefaultFigurePosition', [1   705   560   420]);
   initialize;
   bfinitialize;
   
   verbose = true;
end


%% Read Data
if false
   %%
   
   dirname = '/home/ckirst/Media/ChristophsData/Science/Projects/StemCells/Experiment/Other/Aryeh/Citrine/p1';
   [fileformat, range] = tagformat(dirname, {'z', 't'});
   fileformat
   
   zrange = unique(range(1,:));
   trange = unique(range(2,:));
   
   [min(zrange), max(zrange)]
   [min(trange), max(trange)]
   
   minz = min(zrange);

   %%
   filename = tags2name(fileformat, [3, 1]);
 
   xrange = 500:1012;
   yrange = 500:1012;
   
   t = 40;
   
   img = zeros(length(xrange), length(yrange), length(zrange));
   for z = zrange
      filename = tags2name(fileformat, [z, t]);
      img(:,:,z-minz+1) = imread(fullfile(dirname, filename), 'PixelRegion',{ [500, 1012], [500, 1012]});
   end
  
   img = mat2gray(double(img));
   
   size(img)   
   
   %
   %figure(1); clf
   %implot(img)
   
   %
   figure(1); clf
   set(gcf, 'Name', ['Raw Stack: ' filename ' channel: 1']);  
   downsamplexy = 1;
   imgres = double(img(1:downsamplexy:end, 1:downsamplexy:end, 1:1:end));
   implot3d(imgres)
   
   
   figure(2); clf
   set(gcf, 'Name', ['Raw Stack: ' filename ' channel: 1']);  
   implottiling(imgres)  
else
   filename = '';
end


%% Thresholding / Masking 

th = 2^-2.25% based on local max intensitites of raw image

imgth = img;
imgth(imgth < th) = 0;

imgmask = imgth > 0;
imgmask = imopen(imgmask, strel('disk', 2));


if verbose
   %%
   figure(3)
   clf
   set(gcf, 'Name', ['Thresholded Stack: ' filename ' channel: 1']);
   implottiling(imgth);
   
   fprintf('Masking: foreground: %f / 1.0\n', total(imgmask) / numel(imgmask))
   
   figure(3)
   set(gcf, 'Name', ['Mask: ' filename])
   implottiling(imoverlay(img, imgmask,'red', true));
   
end


%% Filtering

% gaussian smoothing
%imgf = gaussianFilter(imgd,3,10);

% mean shift 
%imgf = meanShiftFilter(imgd, 3, 0.1);

%median
imgf = mat2gray(medianFilter(img, [3, 3, 3]));


imgfc = imgf;
imgfc(imgf > 0.6) = 0.6;
imgfc = mat2gray(imgfc);

if verbose
   %%
   figure(4); clf
   colormap(gray)
   set(gcf, 'Name', ['Filtered Stack: ' filename ' channel: 1']);
   implottiling(imgf);
   
   figure(5); clf
   colormap(gray)
   set(gcf, 'Name', ['Filtered Cutoff Stack: ' filename ' channel: 1']);
   implottiling(imgfc);
   
   
   %figure(5)
   %hist(imgmed(:), 256)
end

%% Prepare Seeding

imgs = mat2gray(imgf);

% add gradient
zsize = size(imgs,3);
%imgg = zeros(size(img));
for z = zsize:-1:1
   %imgg(:,:,z) = imgradient(imgs(:,:,z));
%   imgs(:,:,z) = imgs(:,:,z) - 0.0 * mat2gray(imgradient(imgs(:,:,z)));
end
%imgg = mat2gray(imgg);
%imgs = mat2gray(imgs);


if false
   %%
   figure(60); clf
   colormap(gray)
   set(gcf, 'Name', ['Gradient: ' filename ' channel: 1']);
   implottiling(imgg);
end


% shape based filtering
%imgsf = logFilter(max(imgs(:)) - imgs, [10,10,10], [], 0);
%imgsf = dogFilter(imgs, [12, 12, 4], [] ,[], 0);
%imgsf = diskFilter(imgs, [20, 20, 4], 1, 1, -1); % disk: [ksize or outer-box, width of ring, wt-inner disk, wt-ring]
imgsf = sphereFilter(imgs, [7,7,5], [], 0);

%imgsf = diskFilter(mat2gray(imgg), [16, 16, 7], 7, -0.5, 1);

% post processing
%imgsf(imgsf < 0 ) = 0;
%imgsf = imgsf - min(imgsf(:));
imgsf = mat2gray(imgsf);


if verbose
   figure(6); clf
   colormap(gray)
   set(gcf, 'Name', ['Preprocessed Seeding Image: ' filename ' channel: 1']);
   implottiling(imgsf);
end


% Seeing

%imgmax = imregionalmax(imgf);
imgmax = imextendedmax(imgsf, 0.01);
imgmax = imdilate(imgmax, strel('disk',1));

if verbose
   %%
   figure(10); clf
   set(gcf, 'Name', ['Seeding / Filtered: ' filename ' channel: 1']);
   implottiling(imoverlaylabel(imgsf, bwlabeln(imgmax)));
   
   figure(11); clf
   set(gcf, 'Name', ['Seeding / Raw: ' filename ' channel: 1']);
   implottiling(imoverlaylabel(imgf, bwlabeln(imgmax)));
   

   
end


%% Watershed Segmentation on Image

%dilating the maxima can help
%imgmaxd = imdilate(imgmax, strel('disk', 1));
imgmaxd = imfill(imgmax, 'holes');
%imgmaxd = imgmax;

% watershed
imgmin = imimposemin(max(imgf(:)) - imgf, imgmaxd);
imgws = watershed(imgmin);
imgseg = immask(imgws, imgmask);

%propagation
%imgseg = segmentByPropagation(imgmedf, bwlabeln(imgmaxd), imgmask);

if verbose
   %%
   figure(12); clf
   set(gcf, 'Name', ['Seeding: ' filename ' channel: 1']);
   implottiling(imoverlaylabel(imgf, imgseg, true));
  
end



%% Postprocess Segmentation
% clean up segmentation and alternative diagnositcs
% 

param = setParameter('volume.min',    50,...     % minimal volume to keep (0)
                     'volume.max',    inf,...    % maximal volume to keep (Inf)
                     'intensity.min', -inf, ...  % minimal mean intensity to keep (-Inf)
                     'intensity.max', inf, ...   % maximal mean intensity to keep (Inf)
                     'boundaries',    false, ... % clear objects on x,y boundaries (false)
                     'fillholes',     true,...   % fill holes in each z slice after processing segments (true)
                     'relabel',       true);     % relabel from 1:nlabelnew (true)

[imgpost, stats] = postProcessSegments(imgseg, param);

if verbose
   %%
   figure(100 + figure_offset); clf
   set(gcf, 'Name', ['Final Overlay: ' filename ' channel: 1']);
   implottiling(imoverlaylabel(imgf, imgpost, true))

   figure(101 + figure_offset); clf
   set(gcf, 'Name', ['Final Outline: ' filename ' channel: 1']);
   implotlabeloutline(img, imgpost);
end



end