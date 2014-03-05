% run this file cell by cell to experiment with segmentation parameters. When satisfied
% comment out the display features and unwanted routines... To help others,
% include in ./TestImages a small sample image with a unique and
% descriptive name.
%% preprocessing

img = imread('~/Desktop/martyn/matlab/test_z11.tiff');
imgd = double(img);
imgd = mat2gray(imgd);

imglogvals = log2(imgd(:)+eps);
% ie lower limit 1 on 16 bit image.
imglogvals(imglogvals < -15) = -15;
imglogvals(imglogvals > 0) = 0;


figure(1)
clf
imsubplot(1,2,1);
imshow(imgd);

subplot(1,2,2);
hist(imglogvals, 256);

% all 3 img* files passed to subsequent blocks, nothing to adjust here.

%% thresholding / masking 

param.filter.median = 5; %3;
imgmed = medianFilter(imgd, param.filter.median);

% param.filter.meanshift = 0.1;
% imgmed = meanShiftFilter(imgd, 3, param.filter.meanshift);

thlog = 2^thresholdMixtureOfGaussians(imglogvals, 0.5);
thmed = thresholdMixtureOfGaussians(imgmed, 0.5);

imgth = imgd;
imgth(imgth < thmed) = 0;

imgmask = imgth > 0;
% remove small fragments and small isolated minima with imopen
imgmask = imopen(imgmask, strel('disk', 2));

fprintf('threshold: log(img)= %7.5f, med filter= %7.5f (where max img==1)\n', thlog, thmed);

figure(2)
clf
imsubplot(1,2,1);
imshow(imgth);
imsubplot(1,2,2)
imshow(imoverlay(imgd, imgmask))

% plot statistics to find threshold manually if Mixture of Gaussians fails

figure(3)
subplot(2,3,1)
hist(imgd(:), 256);
title('intensity')

subplot(2,3,2)
hist(imglogvals, 256);
title('log2 intensity')

subplot(2,3,3)
hist(log2(imgmed(:)+eps), 256);
title('log2 median filtered')

% x-axis is effectively number of pixels <= ordinate in following plots
subplot(2,3,4)
plot(sort(imgd(:)))

subplot(2,3,5)
plot(sort(imglogvals))

subplot(2,3,6)
plot(sort(log2(imgmed(:)+eps)))

% imgmask, imgmed passed to subsequent blocks, get mask correct here. 

%% Filtering / Seeding

% gaussian smoothing
%imgdg = gaussianFilter(imgd,3,10);
imgdg = img;

% median filter  ?? also computed above
% imgf = medianFilter(imgdg, 3);

% mean shift 
%imgf = meanShiftFilter(imgd, 3, 0.1);

% disk [ksize or outer-box, width of ring, wt-inner disk, wt-ring]
imgf = diskFilter(imgmed, 7, 1, 1, -1);

% Laplacian of Gaussians
param.filter.logsize = [8, 8];
%imgf = logFilter(max(imgdg(:)) - imgdg, param.filter.logsize);


% find maxima using h-max detection 
param.filter.hmax = 0.02;  %0.02;
%imgmax = imregionalmax(imgf);
imgmax = imextendedmax(mat2gray(imgf), param.filter.hmax);

% Combine nearby points and shrink to single points
% imgmax = imdilate(imgmax, strel('disk', 3));
% imgmax = bwmorph(imgmax,'shrink',inf);


% plot the results. Ideally have one seed per nucleus
figure(12)
imshow([imoverlay(imgd, imgmax) imoverlay(mat2gray(imgf), imgmax)])
title('seeds for watershed overlayed on img and filtered img');

% imgmax passed to segmentation blocks.

%% Watershed Segmentation on Image

%dilating the maxima can help
imgmaxd = imdilate(imgmax, strel('disk', 1));
%imgmaxd = imgmax;

% watershed
imgmin = imimposemin(max(imgmed(:)) - imgmed, imgmaxd);
imgws = watershed(imgmin);
imgseg = double(imgws) .* double(imgmask);

% watershed on gradient
imgmaxd = imdilate(imgmax, strel('disk', 3));
imggrad = imgradient(gaussianFilter(imgd, 3, 10) .* imgmask);
rm = imimposemin(imggrad, imgmaxd);
ws = watershed(rm);
imgsegg = double(ws) .* double(imgmask);


figure(20)
clf
imsubplot(1,2,1)
imshow(imcolorize(imgseg))
title('watershed on image')

imsubplot(1,2,2)
imshow(imcolorize(imgsegg));
title('watershed image gradient')


figure(21)
imshow(immultiply(double(imcolorize(imgseg)), gray2rgb(imgd)))
title('watershed on img overlaid on img')


%% Watershed Segmentation on Image + Gradient

imgmedgrad = imgradient(imgmed);
mi = max(imgmed(:));
mg = max(imgmedgrad(:));
imgmix = imsubtract(imgmed, 1. * (mi/mg) * imgmedgrad);

imgmin = imimposemin(max(imgmix(:)) - imgmix, imdilate(imgmax, strel('disk',3)));
imgws = watershed(imgmin);
imgseg = double(imgws) .* double(imgmask);

figure(30)

ax(1) = imsubplot(2,3,1);
imshow(imoverlay(imgmed, imgmax))

ax(2) = imsubplot(2,3,2);
imshow(imgmedgrad)

ax(3) = imsubplot(2,3,3);
imshow(imgmix);

ax(4) = imsubplot(2,3,4);
imshow(imgmin)

ax(5) = imsubplot(2,3,5);
imshow(imcolorize(imgws));

ax(6) = imsubplot(2,3,6);
imshow(imoverlay(double(imcolorize(imgws)) .* gray2rgb(imgd), imgmax));
linkaxes(ax, 'xy');


%% Segmentation by Propagation

% mixture of image / gradient can improve result
%imgmedgrad = imgradient(imgmed);
%mi = max(imgmed(:));
%mg = max(imgmedgrad(:));
%imgprop = imadd(imgmed, 5.0 * mi / mg * imgmedgrad);
imgprop = imgmed;

imgmaxlabel = bwlabel(imgmax);

param.propagation.lambda = 0.2; % weight between pixel intensity (lambda = 0) changes and spatial distance (lambda = 1)
param.propagation.ksize = 1;    % box width for calculating intensity differences
param.propagation.cutoff.distance = Inf;

[imgproplabel, dist] = segmentByPropagation(imgprop, imgmaxlabel, imgmask, param.propagation.lambda, param.propagation.ksize);
imgproplabel(dist>param.propagation.cutoff.distance) = 0;

figure(40)
ax(1) = imsubplot(1,3,1);
imshow(imoverlay(imgprop, imgmax))
ax(2) = imsubplot(1,3,2);
%imshow(imcolorize(imgproplabel))
imshow(double(imcolorize(imgproplabel)) .* gray2rgb(imgd))

distinf0 = dist;
distinf0(dist==Inf) = 0;
ax(3) = imsubplot(1,3,3);
%imshow(imgmask)
imshow(mat2gray(distinf0))

% trick to make all axes scale together.
linkaxes(ax, 'xy')

%figure(41)
%imshow([imgmedgrad, imgprop])

%% Ray Segmentation - see Test/testRaySegmenetation.m

edit ./Test/testRaySegmentation.m

%% Clean up segmentation and alternative diagnositcs

lbl = bwlabel(imgseg > 0);
stats = regionprops(logical(lbl), 'Area', 'PixelIdxList');
min_area = 40;
keep = [stats.Area] >= min_area;
% remove small nucs from lbl
for i = find(~keep)
    lbl(stats(i).PixelIdxList) = 0;
end

[bndry, lbl] = bwboundaries(lbl > 0, 8, 'noholes');
stats = regionprops(logical(lbl), 'Area', 'Centroid', 'PixelIdxList');
imgmask = lbl > 0;

edges = false(size(imgseg));
for i = 1:length(bndry)
    edges( sub2ind(size(imgseg), bndry{i}(:,1), bndry{i}(:,2)) ) = 1;
end
fprintf('Identified %d nuc in segmentation after discarding %d < min-area Min,Max area= %d %d\n',...
    length(stats), sum(~keep), min([stats.Area]), max([stats.Area]) );
figure(51), imshow(imoverlay(imgmed, edges, [1,0,0]) )

%% ROI: Thresholding and Closing / Opening Original Image - Tests

threshold = 3000/ 65535.;
imgroi = im2bw(imgmed, threshold);


% larger sizes exclude dividing cells !!
imgroi = imerode(imgroi, strel('disk',5));
imgroi = imdilate(imgroi, strel('disk',5));
%imgroi = bwmorph(imgroi, 'open', 15);
%imgroi = bwmorph(imgroi, 'close');

imgthres = zeros(size(imgd));
imgthres(imgroi) = imgd(imgroi);

imgtmed = medianFilter(imgthres,3);
param.filter.median = 5; %3;
imgmed = medianFilter(imgd, param.filter.median);

% param.filter.meanshift = 0.1;
% imgmed = meanShiftFilter(imgd, 3, param.filter.meanshift);

thlog = 2^thresholdMixtureOfGaussians(imglogvals, 0.5);
thmed = thresholdMixtureOfGaussians(imgmed, 0.5);

imgth = imgd;
imgth(imgth < thmed) = 0;

imgmask = imgth > 0;
% remove small fragments and small isolated minima with imopen
imgmask = imopen(imgmask, strel('disk', 2));

fprintf('threshold: log(img)= %7.5f, med filter= %7.5f (where max img==1)\n', thlog, thmed);

figure(2)
clf
imsubplot(1,2,1);
imshow(imgth);
imsubplot(1,2,2)
imshow(imoverlay(imgd, imgmask))

% plot statistics to find threshold manually if Mixture of Gaussians fails

figure(3)
subplot(2,3,1)
hist(imgd(:), 256);
title('intensity')

subplot(2,3,2)
hist(imglogvals, 256);
title('log2 intensity')

subplot(2,3,3)
hist(log2(imgmed(:)+eps), 256);
title('log2 median filtered')

% x-axis is effectively number of pixels <= ordinate in following plots
subplot(2,3,4)
plot(sort(imgd(:)))

subplot(2,3,5)
plot(sort(imglogvals))

subplot(2,3,6)
plot(sort(log2(imgmed(:)+eps)))

% imgmask, imgmed passed to subsequent blocks, get mask correct here. 
imgtmsf = meanShiftFilter(imgthres, 3, 0.1);


figure(100
subplot(2,1,1)
imshow([imgd, imgthres])
subplot(2,1,2)
imshow([imoverlay(imgd, imgroi), imoverlay(imgthres, imgroi)])



