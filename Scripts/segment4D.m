function [stacksegmented, segmentprops] = segment4D(filename, seriesid)
% script for segmentation 

% Segmentation in 4D

verbose = true;

%% Init

close all
clear all
clc

set(0, 'DefaultFigurePosition', [1   705   560   420]);

initialize()

ijinitialize();

%% Data file

filename = '/home/ckirst/Science/Simulation/Matlab/StemCell3D/Test/Images/Develop/Wnt/wnt_clone8_again_feb09.lif';
seriesid = 1;

%% Load Data

%stackraw = imread_bf();
lifdata = imread_bf(filename, seriesid);
stackraw = lifdata(:,:,:,1,1);

%% Load Single Image Data

%dir(fullfile(datadir, 'W1F188T0001*C1.tif'))

boxr = [1 1 1];

if verbose
   figure(1)
   clf
   set(gcf, 'Name', ['Raw Stack: ' filename ' series:' num2str(seriesid)]);
   implot3d(stackraw, 'BoxRatios', boxr, 'Texture', 'z');
end

%%

%ijplot3d(stackraw, 'PixelDepth', 1)

%% threshold

th = 2^thresholdMixtureOfGaussians(log2(stackraw(stackraw>0)), 0.1);
stackth = stackraw;
stackth(stackth < th) = 0;
stackmask = stackth > 0;
stackmask = imdilate(stackmask, strel('disk', 3));
stackmask = imclose(stackmask, strel('disk', 2));


stackmed = medianFilter(stackraw, [3, 3, 3]);
stackmed = stackmed .* stackmask;


if verbose
   figure(2)
   clf
   set(gcf, 'Name', ['Thresholded: ' filenames]);
   ax(1) = imsubplot(1,2,1);
   implot3d(stackmask, 'BoxRatios', boxr, 'Texture', 'z' );
   ax(2) = imsubplot(1,2,2);
   implot3d(stackmed, 'BoxRatios', boxr, 'Texture', 'z' );
   %linkprop(ax, {'CameraPosition', 'CameraViewAngle', 'XLim', 'YLim', 'ZLim'});
   
   figure(3)
   set(gcf, 'Name', ['Mask: ' filenames]);
   immontage(stackmask);
   
   figure(4)
   set(gcf, 'Name', ['Masked Median: ' filenames]);
   immontage(stackmed);
end

%% Seeding

% improve seeding by adding gradients

stackgrad = zeros(size(stackmed));
for z = 1:size(stackmed,3)
   stackgrad(:,:,z) = imgradient(stackmed(:,:,z));
end

%figure(103)
%clf
%implot3d(stackgrad)

mi = max(stackmed(:));
mg = max(stackgrad(:));
stackmix = stackmed - 0.0 * mi / mg * stackgrad;
stackmix(stackmix < 0) = 0;

if verbose
   
   figure(5)
   clf
   set(gcf, 'Name', ['Mixture: ' filenames]);
   implot3d(stackmix);
   %ijplot3d(stackmix)
   
end

%%

% stacklog = logFilter(iminvert(stackmix), [10,10,10]);
% 
% figure(6)
% set(gcf, 'Name', ['LoGFilter: ' filenames])
% clf
% implot3d(imclip(stacklog - mean(stacklog(:))))


%%
% stackdog = dogFilter(stackmix, [10,10,10]);
% 
% figure(7)
% set(gcf, 'Name', ['DoGFilter: ' filenames])
% clf
% implot3d(imclip(stackdog - mean(stackdog(:))))

%%

% stackdisk = sphereFilter(stackmix, 10);
% 
% figure(8)
% set(gcf, 'Name', ['DiskFilter: ' filenames])
% clf
% implot3d(imclip(stackdisk - mean(stackdisk(:))))


%% select fitlered stack and detect seeds

%stackf = imclip(stackdisk - mean(stackdisk(:)));
%stackf = stackf .* cast(stackmask, class(stackf));
stackf = stackmed;
stackf = diskFilter(stackf, 6);
stackf = log(mat2gray(stackf)+eps) + 5;
stackf(stackf <0) = 0;
stackf = mat2gray(stackf);

stackmax = imextendedmax(mat2gray(stackf), 0.01);
%stackmax = imregionalmax(stackf);
stackmax = stackmax .* cast(stackmask, class(stackmax));

%
%stackmax = imerode(stackmax, strel('disk',3));
%stackmax = imdilate(stackmax, strel('disk',3));


if verbose
   figure(10)
   clf
   set(gcf, 'Name', ['PreSeeding: ' filenames]);
   is = stackf - mean(stackf(:));
   is(is < 0) = 0;
   implot3d(is);
   
   figure(11)
   clf
   set(gcf, 'Name', ['Seeding: ' filenames]);
   implot3d(mat2gray(stackmax));
   
   % overlay in imagej
   stackc = gray2rgb(stackmax);
   %size(stackc);
   stackc(:,:,:,2:3) = 0;
   stackc = stackc + gray2rgb(mat2gray(stackraw));
   
   
   ijplot3d(stackc);
   
end

%% Propagate Seeds

% param.lambda = 0.005;
% param.cutoff.difference = 0.05;
% param.averaging.ksize = 1;
% 
% stackf = stackmed;
% stackf = diskFilter(stackf, 3);
% stackf = mat2gray(stackf);
% %stackf = gaussianFilter(stackf, 7);
% %stackf = medianFilter(stackf, 5);
% 
% %stackf = gaussianFilter(stackf, 5);
% stackf = log(stackf + eps) + 5;
% stackf(stackf < 0) = 0;
% stackf = mat2gray(stackf);
% %stackf = log2(stackf + eps) + 5;
% %stackf(stackf < 0) = 0;
% 
% stacklabel = bwlabeln(stackmax);
% 
% %ijplot3d(imcolorize(stacklabel))
% 
% stackprop = seedPropagation(stackf,stacklabel, stackmask, param);
% 
% 
% figure(66)
% clf
% implot3d(stackprop);
% 
% 
% %stackc = gray2rgb(stackprop > 0);
% %size(stackc);
% %stackc(:,:,:,2:3) = 0;
% stackc = imcolorize(stackprop);
% stackc = stackc + gray2rgb(mat2gray(stackraw));
% 
% 
% ijplot3d(stackc)


%% Segment


% 3d water shed on image + gradient 
mi = max(stackmed(:));
mg = max(stackgrad(:));
stackws = stackmed - 0.25 * mi / mg * stackgrad;
stackws(stackws < 0) = 0;

%stackws = mat2gray(stackmed) + 0 * mat2gray(stackf);
rm = imimposemin(iminvert(stackws), stackmax);

ws= watershed(rm);
ws = ws .* cast(stackmask, class(ws));


if verbose
   %figure(13)
   %set(gcf, 'Name', ['Segmentation: ' filenames])
   %clf
   %implot3d(ws)
   
   
   %figure(14)
   %clf
   %set(gcf, 'Name', ['Segmentation: ' filenames])
   %imsurfaceplot3d(ws)
   
   %ijplot3d(imcolorize(ws) + gray2rgb(mat2gray(stackraw)))
   
end


%% postprocess

seg = double(ws) .* double(stackmask);

seg2 = bwlabeln(seg > 0);

%length(unique(seg(:)))
%length(unique(seg2(:)))


% segments size

area = regionprops(seg2, 'area');
area = [area.Area];

%px = regionprops(seg2, 'PixelIdxList');
%px = {px.PixelIdxList};

idx = find(area < 200);
segth = seg2;
for i = idx
   segth(segth == i) = 0;
end

% fill possible holes 

segh = segth;
for s = 1:size(segth,3)
        segh(:,:,s) = imfill(segh(:,:,s), 'holes');
end

% remove boundary objects

% we dont want to clear objects that border in z
pseg = zeros(size(segh) + [0 0 2]);
pseg(:, :, 2:end-1) = segh;
pseg = imclearborder(pseg);
segh = pseg(:,:,2:end-1);
segh = bwlabeln(segh>0);

props = regionprops(segh, stackraw, 'Area', 'MaxIntensity', 'MinIntensity','MeanIntensity', 'PixelValues', 'Centroid');


if verbose
   % plotting
   %figure(30)
   %clf
   %imshow3d(seg2);
   
   
   %figure(32)
   %clf
   %imshow3d(segh)
   
   ijplot3d(imcolorize(segh) + gray2rgb(mat2gray(stackraw)));
   

   figure(33)
   nbin = 25;
   subplot(2,3,1)
   hist([props.Area], nbin)
   title('Area [Pixel]')
   subplot(2,3,2)
   hist([props.MaxIntensity], nbin)
   title('Max Intensity')
   subplot(2,3,3)
   hist([props.MinIntensity], nbin)
   title('Min Intensity')
   subplot(2,3,4)
   hist([props.MeanIntensity], nbin)
   title('Mean Intensity')
   
end


%% return value

stacksegmented = segh;
segmentprops = props;

end


