%%%%%%%%%%%%%%%%%%%%%
%%% Test Stiching %%%
%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
initialize

%%
addpath('./Stitching/Test')

%%

%%%%%%%%%%%%%%%%%%%%%
%%% 2D - 2 Images %%%
%%%%%%%%%%%%%%%%%%%%%

%% Load and Split a Test Image

img = imread('./Test/Images/hESCells_DAPI.tif');
img = mat2gray(img);
size(img)

%
%horizontal
img1 = img(1:300, 1:end-30);
img2 = img(512-300+1:end, 30+1:end);

img1 = img(1:300, 30+1:end);
img2 = img(512-300+1:end, 1:end-30);

size(img1)
size(img2)


figure(1); clf
implot(img)

figure(2); clf
imsubplot(1,2,1)
implot(img1);
imsubplot(1,2,2);
implot(img2);

alignfigures


%% alignment 

shift = align2ImagesByRMS(img1, img2)

figure(2); clf;
plot2AlignedImages(img1, img2, shift)


%%

shift = align2ImagesByRMS(img2, img1)

figure(2); clf;
plot2AlignedImages(img2, img1, shift)


%%

%
%vertical
img1 = img(:,1:290);
img2 = img(:, 512-290+1:end);

size(img1)
size(img2)


figure(1); clf
implot(img)

figure(2); clf
imsubplot(1,2,1)
implot(img1);
imsubplot(1,2,2);
implot(img2);


%% alignment 

shift = align2ImagesByRMS(img1, img2)

figure(2); clf;
plot2AlignedImages(img1, img2, shift)


%%

shift = align2ImagesByRMS(img2, img1)

figure(2); clf;
plot2AlignedImages(mat2gray(img2), mat2gray(img1), shift)



%% Load two images

img1 = imread('./Test/Images/hESCells_tiling/W1F034T0001Z05C1.tif');
img2 = imread('./Test/Images/hESCells_tiling/W1F035T0001Z05C1.tif');

img1 = impqlpermute(img1, 'yx', 'pq');
img2 = impqlpermute(img2, 'yx', 'pq');

size(img1)

figure(1); clf
imsubplot(1,2,1);
implot(double(img1));
imsubplot(1,2,2);
implot(img2);





%% alignment 

shift = align2ImagesByRMS(img1, img2, setParameter('overlap.min', 50))

figure(2); clf;
plot2AlignedImages(mat2gray(img1), mat2gray(img2), shift)


%%

shift = align2ImagesByRMS(img2, img1, setParameter('overlap.min', 70))

figure(2); clf;
plot2AlignedImages(mat2gray(img2), mat2gray(img1), shift)




%% Test algin all images



%% Load and Split a Test Image

img = imread('./Test/Images/hESCells_DAPI.tif');
img = mat2gray(img);
size(img)

%
%horizontal
img1 = img(1:300, 1:300);
img2 = img(512-300+1:end, 1:300);
img3 = img(1:300, 512-300+1:end);
img4 = img(512-300+1:end, 512-300+1:end);

var2char({size(img1), size(img2), size(img3), size(img4)})


figure(1); clf
implot(img)

figure(2); clf
imsubplot(2,2,1)
implot(img3);
imsubplot(2,2,2);
implot(img4);
imsubplot(2,2,3)
implot(img1);
imsubplot(2,2,4);
implot(img2);


%%
sh = alignImages({img3, img4; img1, img2});

var2char(sh)

figure(2); clf;
plotAlignedImages({img3, img4; img1, img2}, sh)

