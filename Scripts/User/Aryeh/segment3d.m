function imgseg = segment3d(img, verbose)
% 3d segmentation for Wnt data

%% Init

if false
   %%
   close all
   clear all
   clc
   set(0, 'DefaultFigurePosition', [1   705   560   420]);
   
   initialize()
   
   %ijinitialize();
   
   %%
   imarisinitialize();
end

if nargin < 2
   %%
   verbose = true;
   
   %%
   verbose = false;
end

%% get some test data

if false
   %%
   filename = '/home/ckirst/Science/Simulation/Matlab/StemCell3D/Test/Images/Develop/Aryeh/140305_RUES2_36hBMP4_Bra_Snail_Sox2.lif';
   %lifdata = imread_bf(filename, struct('series', 21, 'channel', 1));
  
   lifdata = imread_bf(filename, struct('series', 21, 'time', 1, 'channel', 1, 'y', 1000 + [0, 512]));
   %lifdata = imread_bf(filename, struct('series', 2, 'time', 1, 'channel', [], 'x', [1, 126], 'y', [1,126]));
   lifdata = imzinvert(lifdata);
   
   img = lifdata(:,:,:,1,1);
   
   clear lifdata
else
   filename = '';
end

%% Data from imaris

if false
    %% 
    imarisstart();
    
    %%
    img = imarisget('Volume', 0,0); 
end


%%
if verbose
   figure(8)
   clf
   imgres = img(1:5:end, 1:5:end, 1:2:end);
   implot3d(imgres)
end



%% initialize / prefiltering

imgd = double(img);
imgd = mat2gray(imgd);

%imglogvals = log2(imgd(:)+eps);
%imglogvals(imglogvals < -5) = -5;
%imglogvals(imglogvals > 0) = 0;

if verbose
   figure(1)
   clf
   set(gcf, 'Name', ['Raw Stack: ' filename ' channel: 1']);
   implot3d(imgd);
   
   %figure(2)
   %subplot(1,2,1);
   %hist(imglogvals, 256)
   %subplot(1,2,2);
   %hist(imgd(:), 256);
end


%% filter 
param.filter.median.ksize = [3, 3, 3]; %3;
imgmed = mat2gray(medianFilter(imgd, param.filter.median.ksize));
%imgmed = meanShiftFilter(imgmed, 4, 0.1);


imgmedf = imgmed;
imgmedf(imgmed > 0.6) = 0.6;
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
%th = thresholdMixtureOfGaussians(imgmedf, 0.5)
%th = thresholdEntropy(imgmed);
%th = thresholdMutualEntropy(imgmed);
%th = thresholdOtsu(imgmed);
th = 0.175;

imgth = imgmedf;
imgth(imgth < th) = 0;

imgmask = imgth > 0;
% remove small fragments and small isolated minima with imopen
imgmask = imopen(imgmask, strel('disk', 3));


if verbose
   %%
   figure(4)
   clf
   set(gcf, 'Name', ['Thresholded Stack: ' filename ' channel: 1']);
   implot3d(imgth);
   
   %figure(5)
   %clf
   %set(gcf, 'Name', ['Mask: ' filename ' channel: 1']);
   %imsurfaceplot3d(imgmask);
   
   %ijplot3d(imgth, 'PixelDepth', 5)
   
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
imgf = mat2gray(imgf);


if verbose
   figure(17)
   clf
   imsubplot(1,2,1)
   implot3d(imgf)
end

%param.filter.logsize = [10,10,10];
%imgf = logFilter(max(imgf(:)) - imgf, param.filter.logsize, [], 0);

imgf = dogFilter(imgf, [15, 15, 7], [] ,[], 0);


if verbose 
   imsubplot(1,2,2)
   imgs = imgf;
   imgs(imgs< 0) = 0;
   implot3d(mat2gray(imgs))
end


%imgf = sphereFilter(imgf, [7,7,7], [], 0);

%imgf = diskFilter(imgf, [11,11,9], 1, 1, -1, 0);
%imgf(imgf < 0 ) = 0;
imgf = imgf - min(imgf(:));
imgf = mat2gray(imgf);




%%%
% find maxima using h-max detection 
param.filter.hmax = 0.01;  %0.02;
%imgmax = imregionalmax(imgf);
imgmax = imextendedmax(mat2gray(imgf), param.filter.hmax);
imgmax = imdilate(imgmax, strel('disk',2));

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
%imgmaxd = imdilate(imgmax, strel('disk', 1));
imgmaxd = imfill(imgmax, 'holes');
%imgmaxd = imgmax;

% watershed
imgmin = imimposemin(max(imgmedf(:)) - imgmedf, imgmaxd);
imgws = watershed(imgmin);
imgseg = immask(imgws, imgmask);

%propagation
%imgseg = segmentByPropagation(imgmedf, bwlabeln(imgmaxd), imgmask);



% %% Clean up segmentation and alternative diagnositcs
% 
imgseg = imlabelseparate(imgseg);
stats = regionprops(imgseg, 'Area', 'PixelIdxList');

min_area = 40;
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
   %implot3d(imgseg);
   imsurfaceplot3d(imgseg)
end


%%
if verbose
   %%
   figure(23)
   clf
   implotlabeloutline(img, imgseg)
end


%%
if true
   ijplot3d(imcolorize(imgseg) + gray2rgb(mat2gray(img)), 'PixelDepth', 5);
end


%%

size(imlabel(imgseg))

%% calculate surfaces and move to imaris

[surf, fac, norm] = imsurface(imgseg, 'all');

%%
sf = {surf, fac, norm};


%% save this stuff

%save('./Test/Images/Develop/Aryeh/140305_RUES2_36hBMP4_Bra_Snail_Sox2_segmetation_imaris.mat', 'imgseg')
%save('./Test/Images/Develop/Aryeh/140305_RUES2_36hBMP4_Bra_Snail_Sox2_surfaces_imaris.mat', 'sf')


%% push stufff to imaris

nset = size(surf);
sfset = surf(1:nset);
fcset = fac(1:nset);
nmset = norm(1:nset);

imarissetsurface('Aryeh Segments', sfset, fcset, nmset);




%% 

%% save the result

%save('./Test/Images/Develop/Aryeh/140305_RUES2_36hBMP4_Bra_Snail_Sox2_segmetation_2.mat', 'imgseg')

%load('./Test/Images/Develop/Aryeh/140305_RUES2_36hBMP4_Bra_Snail_Sox2_segmetation.mat')



%% 
stats = statisticsSegments(double(img), imgseg);

%%
save('./Test/Images/Develop/Aryeh/140305_RUES2_36hBMP4_Bra_Snail_Sox2_segmetation_2_stats_ch1.mat', 'stats')

%%
 
lifdata = imread_bf(filename, struct('series', 21, 'time', 1, 'channel', 2, 'y', 1000 + [0, 512]));
%lifdata = imread_bf(filename, struct('series', 2, 'time', 1, 'channel', [], 'x', [1, 126], 'y', [1,126]));
lifdata = imzinvert(lifdata);
lifdata = squeeze(lifdata);

%% 
stats2 = statisticsSegments(lifdata, imgseg);

%%
save('./Test/Images/Develop/Aryeh/140305_RUES2_36hBMP4_Bra_Snail_Sox2_segmetation_2_stats_ch2.mat', 'stats2')


%%
 
lifdata = imread_bf(filename, struct('series', 21, 'time', 1, 'channel', 3, 'y', 1000 + [0, 512]));
%lifdata = imread_bf(filename, struct('series', 2, 'time', 1, 'channel', [], 'x', [1, 126], 'y', [1,126]));
lifdata = imzinvert(lifdata);
lifdata = squeeze(lifdata);

%% 
stats3 = statisticsSegments(lifdata, imgseg);

%%
save('./Test/Images/Develop/Aryeh/140305_RUES2_36hBMP4_Bra_Snail_Sox2_segmetation_2_stats_ch3.mat', 'stats3')

%%
 
lifdata = imread_bf(filename, struct('series', 21, 'time', 1, 'channel', 4, 'y', 1000 + [0, 512]));
%lifdata = imread_bf(filename, struct('series', 2, 'time', 1, 'channel', [], 'x', [1, 126], 'y', [1,126]));
lifdata = imzinvert(lifdata);
lifdata = squeeze(lifdata);

%% 
stats4 = statisticsSegments(lifdata, imgseg);

%%
save('./Test/Images/Develop/Aryeh/140305_RUES2_36hBMP4_Bra_Snail_Sox2_segmetation_2_stats_ch4.mat', 'stats4')



%%
if verbose
   %%
   figure(11)
   clf
   %implot3d(imgseg);
   %cdata = [stats.Centroid];
   %cdata = cdata(3,:);
   
   %cdata = [stats.MeanIntensity];
   cdata = [stats.MaxIntensity];
   cdata = [stats.Volume];
   
  
   param.color.data = cdata;
   imsurfaceplot3d(imgseg(1:512, 1:512, :), param)
end

%%
if verbose
   
   xmeasure = 'Volume';
   
   figure(45)
   clf
   
   ymeasures = {'MeanIntensity', 'MinIntensity','MaxIntensity', 'MedianIntensity'};
   for i = 1:4
      ax(i) = subplot(2,2,i);
      plot([stats.(xmeasure)], [stats.(ymeasures{i})], '*');
      xlabel(xmeasure); ylabel(ymeasures{i});
   end
  
   linkaxes(ax, 'x')
end


%%
if verbose
   
   xmeasure = 'z';
   xdata = [stats.Centroid];
   xdata = xdata(3,:);
   
   figure(45)
   clf
   ymeasures = {'MeanIntensity', 'MinIntensity','MaxIntensity', 'MedianIntensity'};
   for i = 1:4
      ax(i) = subplot(2,2,i);
      plot(xdata, [stats.(ymeasures{i})], '*');
      xlabel(xmeasure); ylabel(ymeasures{i});
   end
  
   linkaxes(ax, 'x')
   
end



end

