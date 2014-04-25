%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Example: Triple Reporter - High Resolution Stack %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data: Albert Ruzo

clear all
close all
clc

initialize
ijinitialize;

%% Read Images from Movie
img = imread_avi('./Test/Images/hESC_tripple_stack.avi', 'time', [1 70]);
size(img)

ijplot3d(img, 'PixelDepth', 3)

%% Read  Background from Movie
imgback = imread_avi('./Test/Images/hESC_tripple_stack.avi', 'time', [81 100]);
size(imgback)

ijplot3d(imgback, 'PixelDepth', 3)

%% Background 
%imgbacki = sqrt(sum(double(imgback).^2,4));
imgbacki =sum(double(imgback),4);

figure(1); clf
set(gcf, 'Name', 'Histogram Background');
hist(imgbacki(:), 256)

imgbackm = mean(imgbacki, 3);
imgbackmg = gaussianFilter(imgbackm, [50,50]);
%imgbackmg = imopen(imgbackm, strel('disk', 50));

figure(2); clf
set(gcf, 'Name', 'Background / Guassian filtered Background');
implottiling({imgbackm, imgbackmg})
{min(imgbackm(:)), max(imgbackm(:))}
{min(imgbackmg(:)), max(imgbackmg(:))}


%% Background Correction
%imgi = sqrt(sum(double(img).^2,4));
imgi = sum(double(img),4);
{min(imgi(:)), max(imgi(:))}

imgi = imgi - repmat(imgbackmg, 1, 1, size(imgi,3),1);
imgi(imgi < 0) = 0;

imgi = medianFilter(imgi, [3,3,3]);

%%
imgi = imtophat(imgi, strel('disk', 15));


figure(1); clf; colormap gray
set(gcf, 'Name', 'Image - Background');
implottiling(imgi(:,:,1:5:end, :))


figure(2); clf
set(gcf, 'Name', 'Image - Background');
subplot(1,2,1)
hist(imgi(:), 256)
subplot(1,2,2)
imglog = log2(imgi+eps);
imglog(imglog < -5) = -5;
hist(imglog(:), 256)


%% Thresholding

th = 0.0;
imgth = mat2gray(imgi);
imgth(imgth < th) = th;
imgth = mat2gray(imgth);

figure(1); clf; colormap gray
set(gcf, 'Name', 'Thresholded');
implottiling(imgth(:,:,1:5:end))


%% Seeding 

%%
zmax = 70;
zmaxfig = 25;


%%
%imgs= sphereFilter(imgi, [20,20,15]);
%imgs = mat2gray(imgs);

imgs2 = sphereFilter(imgi, [5,5,5]);
imgs2 = mat2gray(imgs2);

imgs2 = medianFilter(imgs2, [3,3,3]);

imgs3 = imtophat(imgs2, strel('disk', 15));


%%
figure(1); clf; colormap gray
set(gcf, 'Name', 'Seed Filtering 2');
implottiling(imgs2(:,:,1:1:zmaxfig, :))

figure(2); clf; colormap gray
set(gcf, 'Name', 'Seed Filtering 3');
implottiling(imgs3(:,:,1:1:zmaxfig, :))


%%
imgmax = imextendedmax(imgs3(:,:,1:zmax), 0.05);

figure(2); clf; 
set(gcf, 'Name', 'Seeds');
colormap jet
implottiling(imoverlaylabel(mat2gray(imgi(:,:,1:1:zmax)), bwlabeln(imgmax(:,:,1:zmax))))


%% Masking

%%
figure(1); clf
hist(imgs3(:),256)

%%
imgmask = mat2gray(imgs3) > 0.03;
imgmask = imopen(imgmask, strel('disk', 3));
imgmask = imdilate(imgmask, strel('disk', 1));

figure(4); clf
%implottiling(imoverlay(mat2gray(imgs3(:,:,1:zmaxfig)), imgmask(:,:, 1:zmaxfig),'r', false))
implottiling(imoverlay(mat2gray(imgs3(:,:,1:5:zmax)), imgmask(:,:, 1:5:zmax),'r', false))

%% Watershed
imgmin = immask(imgs3, imgmask);
imgmin = imimposemin(max(imgmin(:)) - imgmin(:,:,1:zmax), imgmax);

%imgmin = imimposemin(max(imgi(:)) - imgi(:,:,1:zmax), imgmax);

imgws = watershed(imgmin);
%%
imgseg = immask(imgws, imgmask(:,:,1:zmax));
imgseg = imrelabel(imgseg);

%%
figure(42); clf
set(gcf, 'Name', 'Segmentation');
implotlabeloutline(imgs3(:,:,1:zmaxfig), imgseg(:,:,1:zmaxfig));

%%
figure(43); clf
set(gcf, 'Name', 'Segmentation');
implotlabeloutline(imgs2, imgseg);


%%
ijplot3d(imcolorize(imgseg), 'PixelDepth', 3)


%% Postprocess

%%
max(imgseg(:))
imgsep = imlabelseparate(imgseg);
max(imgsep(:))


%%
param = setParameter('volume.min',    1000,...    % minimal volume to keep (0)
                     'volume.max',    inf,...    % maximal volume to keep (Inf)
                     'intensity.min', -inf, ...  % minimal mean intensity to keep (-Inf)
                     'intensity.max', inf, ...   % maximal mean intensity to keep (Inf)
                     'boundaries',    false, ... % clear objects on x,y boundaries (false)
                     'fillholes',     true,...   % fill holes in each z slice after processing segments (true)
                     'relabel',       true);     % relabel from 1:nlabelnew (true)

[imgpost, stats] = postProcessSegments(imgsep, param);
max(imgpost(:))


%% 
figure(47); clf
set(gcf, 'Name', 'Segmentation - Postprocessed');
implotlabeloutline(imgi, imgpost);


%%
ijplot3d(imcolorize(imgpost) + gray2rgb(mat2gray(imgi)), 'PixelDepth', 3);


%% Measure Fluorescence

imgb = img(:,:,:,3) - repmat(imgbackmg, 1, 1, zmax);
imgb = imgb - imopen(imgb, strel('disk', 20));

stats = imstatistics(imgpost, {'MedianIntensity'}, imgb);

%%
imgcol = imcolorize(imgpost, setParameter('color.data', [stats.MedianIntensity], 'color.map', 'b'));

ijplot3d(imgcol, 'PixelDepth' , 3)


%%

ijplot3d(imgray2color(imgb, 'b'), 'PixelDepth' , 3)



%%
% save('./Test/albert.mat')

%%
%load('./Test/albert.mat')
%initialize













