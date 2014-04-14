%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Segmentation of 3D Images - Wnt Example %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example for 3D image segmentation
%


%%
close all
clear all
clc

set(0, 'DefaultFigurePosition', [1   705   560   420]);

initialize()

ijinitialize();


%% Experiment setup
exp = Experiment('Name', 'Wnt', 'Date', datestr(now), 'Description', 'Wnt',...
                 'BaseDirectory', './Test/Data/Experiment/Wnt', ...
                 'ImageDirectoryName', '../../../Images/mESCells_tif_wnt',...
                 'ResultDirectoryName', '.', ...
                 'ReadImageCommandFormat', 'imread(''<directory>/<file>'')',...
                 'ImageFileFormat', 'Zstack_63x_stepsize_0.5um_2um_z<z,3>_c<c,3>.tif', ...
                 'ImageFormatNames', {'z', 'c'},...
                 'ImageFormatRange', {1:29, 1:3});

exp.Info()

%% Read Stack

for z = 12:25 % exp.ImageFormatRange{1}
   img{1}(:,:, z-11) = exp.ReadImage([z, 1]);
end
size(img{1})

imgr = imresample(img{1}, [3, 3, 3]);
size(imgr)

figure(2); clf
implot3d(mat2gray(imgr));

%% Read Other Channels

for c = 2:3
   
   for z = 12:25 % exp.ImageFormatRange{1}
      img{c}(:,:, z-11) = exp.ReadImage([z, c]);
   end
   size(img{c})
   
   
   resamp = [4, 4, 2];
   imgr = imresample(img{c}, resamp);
   size(imgr)
   
   figure(2 + c); clf
   implot3d(mat2gray(imgr));

end

%% some statistics

figure(6)
for c= 1:3
   subplot(1,3,c)
   hist(double(img{c}(:)), 256)
end
   
%%

th = 370;

for c=1:3
   img{c}(img{c} < th) = 0;
end

figure(7)
for c = 1:3
   imgr = imresample(img{c}, resamp);
   imsubplot(1,3,c);
   implot3d(mat2gray(imgr));
end

%%

figure(8)
clf
implottiling(mat2gray(img{1}))



%% Prefiltering / Denoising

c = 1;
imgf = img{c};

%% adding images to enhance cells

imgf = img{1} + img{2} + img{3};
back = imopen(imgf, strel('disk', 25));
imgf = imgf - back;
imgf(imgf < 0 ) = 0;


%imgf = imsharpen(imgf);
%imgf = medianFilter(img{c}, 3);










%imgd = double(img);
%imgd = mat2gray(imgd);

%imglogvals = log2(imgd(:)+eps);
%imglogvals(imglogvals < -5) = -5;
%imglogvals(imglogvals > 0) = 0;

if verbose
   figure(1)
   clf
   set(gcf, 'Name', ['Preprocess Stack: channel: 1']);
   implottiling(mat2gray(imgf));
   
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



% %% Clean up segmentation and alternative diagnositcs
% 
imgseg = imlabelseparate(imgseg); % CK impose mask before watershed.
stats = regionprops(imgseg, 'Area', 'PixelIdxList');

min_area = 40;  % its actually volume use length(pixel...list)
keep = [stats.Area] >= min_area;
for i = find(~keep)
   imgseg(stats(i).PixelIdxList) = 0;
end
imgseg = imfill(imgseg, 'holes');
imgseg = imrelabel(imgseg);

% 
% [bndry, lbl] = bwboundaries(lbl > 0, 8, 'noholes');
% stats = regionprops(logical(lbl), 'Area', 'Centroid', 'PixelIdxList');
% imgmask = lbl > 0;

if verbose
   %%
   figure(11)
   clf
   implot3d(imgseg);
   %imsurfaceplot3d(imgseg);
end

%%
if true
   ijplot3d(imcolorize(imgseg) + gray2rgb(mat2gray(img)), 'PixelDepth', 5);
end


end

