%% preprocessing

img = imread('./Test/Images/raw.tif');
imgd = double(img);
imgd = mat2gray(imgd);

imglogvals = log2(imgd(:)+eps);
imglogvals(imglogvals < -15) = -15;
imglogvals(imglogvals > 0) = 0;


figure(1)
clf
imsubplot(1,2,1)
imshow(imgd);

subplot(1,2,2)
hist(imglogvals, 256)


%% thresholding / masking 
th = 2^thresholdMixtureOfGaussians(imglogvals, 0.5)

imgth = imgd;
imgth(imgth < th) = 0;

imgmask = imgth > 0;
imgmask = imopen(imgmask, strel('disk', 1));
imgmask = imclose(imgmask, strel('disk', 1));


figure(2)
clf
imsubplot(1,2,1)
imshow(imgth);
imsubplot(1,2,2)
imshow(imoverlay(imgd, imgmask))

% plot statistics to find threshold manually if Mixture of Gaussians fails

imgmed = medianFilter(imgd, 3);

figure(3)
subplot(2,3,1)
hist(imgd(:), 256)
title('intensity')

subplot(2,3,2)
hist(imglogvals, 256)
title('log2 intensity')

subplot(2,3,3)
hist(log2(imgmed(:)+eps), 256)
title('log2 median filtered')

subplot(2,3,4)
plot(sort(imgd(:)))

subplot(2,3,5)
plot(sort(imglogvals))

subplot(2,3,6)
plot(sort(log2(imgmed(:)+eps)))



%% Filtering / Seeding

% gaussian smoothing
%imgdg = gaussianFilter(imgd,3,10);
imgdg = img;

% median filter
%imgf = medianFilter(imgdg, 3);

% mean shift 
%imgf = meanShiftFilter(imgd, 3, 0.1);

% disk
imgf = diskFilter(imgmed, 3, 4, 1, -1);

% Laplcaina of Gaussians
param.filter.logsize = [8, 8];
%imgf = logFilter(max(imgdg(:)) - imgdg, param.filter.logsize);


% find maxima using h-max detection 
param.filter.hmax = 0.02;
%imgmax = imregionalmax(imgf);
imgmax = imextendedmax(mat2gray(imgf), param.filter.hmax);

%%% Combine nearby points and shrink to single points
%imgmax = imdilate(imgmax, strel('disk', 3));
%imgmax = bwmorph(imgmax,'shrink',inf);


% plot the results
figure(10)
imshow([imoverlay(imgd, imgmax) imoverlay(mat2gray(imgf), imgmax)])


%% Watershed Segmentation

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

imsubplot(1,2,2)
imshow(imcolorize(imgsegg));


figure(21)
imshow(immultiply(double(imcolorize(imgseg)), gray2rgb(imgd)))




%% Watershed on Image + Gradient

imgmedgrad = imgradient(imgmed);
mi = max(imgmed(:));
mg = max(imgmedgrad(:));
imgmix = imsubtract(imgmed, .5 * mi / mg * imgmedgrad);

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

linkaxes(ax, 'xy')

%figure(41)
%imshow([imgmedgrad, imgprop])



%% Ray Segmentation - see Test/testRaySegmenetation.m

edit ./Test/testRaySegmentation.m



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
imgtmsf = meanShiftFilter(imgthres, 3, 0.1);


figure(100)
subplot(2,1,1)
imshow([imgd, imgthres])
subplot(2,1,2)
imshow([imoverlay(imgd, imgroi), imoverlay(imgthres, imgroi)])




