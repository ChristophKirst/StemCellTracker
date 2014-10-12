%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test findRegionOfInfluence %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%% Load some image data

img = imread('./Test/Images/hESCells_tif_tiling/W1F038T0001Z05C1.tif');

img = img(1:128, 1:128);

figure(1); clf
implot(img)

%% Detect Centers

imgf = mat2gray(log(double(img+eps)));
imgf = filterMedian(imgf, 3);
imgf(imgf < 0.075) = 0;
imgmask = imgf > 0;
imgf = filterDisk(imgf, 4);

imgmax = imextendedmax(imgf, 0.015);
imgmax = immask(imgmax, imgmask);

implot(imoverlay(imgf, imgmax))

imglab = bwlabeln(imgmax);

%% findRegionOfInfluence

imgroi = findRegionOfInfluence(imglab, setParameter('distance.max', 10));

figure(2); clf
implottiling({imoverlay(img, imgmax), imoverlaylabel(img, imgroi, true)})

%note: in applications start fromsegmented nuclei !


%% check labeling

d = [];
for i = 1:max(imglab(:))
   idx = find(imglab == i);
   d = [d, any(~(imgroi(idx) == i))];
end
d










%% testing ....

dist = bwdist(imgmax);
max(dist(:))


mask = dist < 4;

dist = mat2gray(dist);
ws = watershed(dist);



wsm = immask(ws, mask);
colormap jet
implottiling({imgmax, dist, imcolorize(ws), imcolorize(wsm)})












