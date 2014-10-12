%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test joinSeedsByRays %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%% Load some image data

img = imread('./Test/Images/hESCells_tif_tiling/W1F038T0001Z05C1.tif');

img = img(128:256, 128:256);

figure(1); clf
implot(img)

%% Detect Centers

imgf = mat2gray(log(double(img+eps)));
imgf = filterMedian(imgf, 3);
imgf(imgf < 0.075) = 0;
imgmask = imgf > 0;
imgf = filterDisk(imgf, 3);

imgmax = imextendedmax(imgf, 0.005);
imgmax = immask(imgmax, imgmask);

implot(imoverlay(imgf, imgmax))

imglab = bwlabeln(imgmax);

%% joinSeedsByRays

param = setParameter('threshold.min',        0.075, ... % if profile comes below this absolute intensity objects are different
                     'threshold.max',        0.9,...    % if profile is always above this threshold, objects are joined
                     'threshold.change'    , 0.3, ...   % maximal rel change in intensitiy above objects are assumed to be different
                     'threshold.gradient'  , 0.4, ...   % maximal absolute gradient change below objects are joined, only if gradient image is supplied
                     'cutoff.distance'     , 10, ...    % maximal distance between labels (= 20)
                     'averaging.ksize'     , 3, ...     % ksize to calculate reference mean intensity (=3)
                     'addline'             , true);     % add a line between joined label (true)

[imgjoin, pairs, joins] = joinSeedsByRays(imglab, imgf, param);


figure(3)
implot(imoverlaylabel(img, imgjoin))

figure(4)
clf
imsubplot(1,2,1)
implot(imoverlay(img, imglab));
for i = 1:size(pairs, 1);
   col = hsv2rgb([i/size(pairs,1), 1, 1]);
   plotSeedPairs(imglab, pairs(i, :), col)
end

imsubplot(1,2,2)
implot(imoverlay(img, imglab));
for i = 1:size(joins,1)
   %col = hsv2rgb([i/size(joins,1), 1, 1])
   plotSeedPairs(imglab, joins(i,:),  col)
end













%% testing / snippets....

cent  = imstatistics(imglab, 'Centroid');
cent = round([cent.Centroid]);

isize = size(imglab);

clabel = zeros(isize);
clabel(imsub2ind(isize, cent')) = 1;
%clabel = clabel * imglab;

figure(2)
clf
implot(imoverlay(img, clabel))

%plotSeedPairs(imglab, pairs)


%%


imwrite(mat2gray(double(img)), '~/Desktop/test.jpg')
