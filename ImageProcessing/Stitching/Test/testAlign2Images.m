%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test Image Alignment of 2 Images %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

initialize

%% Load Test Image

img = loadTestImage();
img = mat2gray(img);

size(img)

%% Images with large overlap
img1 = img(1:453, 1:500);
img2 = img(20:end, 58:end);
imgs = {img1; img2};

%img1 = img(1:end-20,:);
%img2 = img(30:end, :);
%imgs = {img1, img2};


figure(1); clf;
implottiling(imgs, 'link', false);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Align 2 Images

%% Optimization

tic
sh = align2ImagesByOptimization(img1, img2)
toc

figure(2);
plot2AlignedImages(img1, img2, sh)


%% RMS of differences

tic
sh = align2ImagesByRMS(img1, img2, 'overlap.min', 1)
toc

figure(2); clf
plot2AlignedImages(img1, img2, sh)


%% Correlation

tic
sh = align2ImagesByCorrelation(img1, img2, 'overlap.min', 1, 'meancorrection', true)
toc

figure(2); clf
plot2AlignedImages(img1, img2, sh)


%% Hugin feature detection

tic
sh = align2ImagesByHugin(img1, img2)
toc

figure(2); clf
plot2AlignedImages(img1, img2, sh)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Aligment 2 Images Left Right


%% Tiles with Small Overlap

img1 = img(1:230,1:end-30);
img2 = img(200:end,30:end);
imgs = {img1; img2};

figure(1); clf;
implottiling(imgs, 'link', false);


%% RMS

tic
sh = align2ImagesLeftRightByRMS(img1, img2, 'overlap.max', 70, 'overlap.min', 10, 'shift.max', 50)
toc 

figure(2); clf
plot2AlignedImages(img1, img2, sh)


%% Correlation

tic
sh = align2ImagesLeftRightByCorrelation(img1, img2, 'overlap.max', 70, 'overlap.min', 10, 'shift.max', 50)
toc 

figure(2);
plot2AlignedImages(img1, img2, sh)


%% Hugin

tic
sh = align2ImagesLeftRightByHugin(img1, img2, 'overlap.max', 80, 'overlap.min', 0, 'shift.max', 50)
toc 

figure(2); clf
plot2AlignedImages(img1, img2, sh)


%% Optimization 

tic
sh = align2ImagesLeftRightByOptimization(img1, img2, 'overlap.max', 80, 'overlap.min', 0, 'shift.max', 150)
toc 

figure(2);
plot2AlignedImages(img1, img2, sh)

% (note: fails here, matalb routine slow and non-robust)












