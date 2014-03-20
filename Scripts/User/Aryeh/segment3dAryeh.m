%function imgseg = segment3d(img, verbose)
%
% 3d segmentation for geometrically confined stemm cell colonies
% 

%% Initialize Segmenter
if false
   %%
   close all
   clear all
   clc
   set(0, 'DefaultFigurePosition', [1   705   560   420]);
end
initialize

%% Init

if nargin < 2
   verbose = false;
end

useimaris = false;
savefile = '';
verbose = false;

%%
if false
   %% Init for Testing
   bfinitialize();   
   initialize();

   %%
   if useimaris
      imarisinitialize();
   end

   %% Data form disk

   filename = '/home/ckirst/Science/Simulation/Matlab/StemCell3D/Test/Images/Develop/Aryeh/140305_RUES2_36hBMP4_Bra_Snail_Sox2.lif';
   
   %lifdata = imread_bf(filename, struct('series', 21, 'channel', 1));
   %lifdata = imread_bf(filename, struct('series', 21, 'time', 1, 'channel', 1, 'y', 1000 + [0, 512]));
   %img = lifdata(:,:,:,1,1);
   %clear lifdata
   
   xr = [1, 512];  % use [] for all
   yr = [1, 512];
   cr = 1;
   se = 2;
   ti = 1;
   img = imread_bf(filename, struct('series', se, 'time', ti, 'channel', cr, 'x', xr, 'y', yr));
   img = imzreverse(squeeze(img));
   
   
   
   %%  figure aryeh
   xr = [1000, 1300];  % use [] for all
   yr = [1300, 1600];
   cr = 1;
   se = 21;
   ti = 1;
   img = imread_bf(filename, struct('series', se, 'time', ti, 'channel', cr, 'x', xr, 'y', yr));
   img = imzreverse(squeeze(img));

   
   %% Data from imaris

   if useimaris
      %%
      imarisstart();
      
      %%
      img = imarisget('Volume', 0,0); 
   end
   
   
else
 
   filename = '';
end


%% Plot low resolution 
if verbose
   %%
   figure(1)
   clf
   downsamplexy = 5;
   imgres = img(1:downsamplexy:end, 1:downsamplexy:end, 1:1:end);
   implot3d(mat2gray(imgres))
end



%% initialize / prefiltering

imgd = double(img);
imgd = mat2gray(imgd);

%imglogvals = log2(imgd(:)+eps);
%imglogvals(imglogvals < -5) = -5;
%imglogvals(imglogvals > 0) = 0;

if verbose
   figure(2)
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
   %imsurfaceplot3d(imgmask);   lnew = lnew + 1;
   
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
   imsubplot(1,2,1);
   implot3d(imgf);
end

%param.filter.logsize = [10,10,10];
%imgf = logFilter(max(imgf(:)) - imgf, param.filter.logsize, [], 0);

imgf = dogFilter(imgf, [15, 15, 7], [] ,[], 0);


if verbose 
   imsubplot(1,2,2);
   imgs = imgf;
   imgs(imgs< 0) = 0;
   implot3d(mat2gray(imgs));
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
   imsurfaceplot3d(imgseg);
end


%%
if verbose
   %%
   figure(23)
   clf
   implotlabeloutline(img, imgseg);
end


%%
if true
   ijplot3d(imcolorize(imgseg) + gray2rgb(mat2gray(img)), 'PixelDepth', 5);
end

%%
size(imlabel(imgseg))
imgseg = imrelabel(imgseg);
size(imlabel(imgseg))

%% Calcualte Statistics 

stats = cell(4,1);
stats{1} = statisticsSegments(double(img), imgseg);

%%
for ch = 2:4

   img2 = imread_bf(filename, struct('series', se, 'time', ti, 'channel', ch, 'x', xr, 'y', yr));
   img2 = imzreverse(squeeze(img2));
   
   stats{ch} = statisticsSegments(img2, imgseg);
end




%% Calculate Surfaces

[surf, fac, norm] = imsurface(imgseg, 'all');
surfaces = {surf, fac, norm};

%% save this stuff

if ~isempty(savefile)
   %%
   %save('./Test/Images/Develop/Aryeh/140305_RUES2_36hBMP4_Bra_Snail_Sox2_segmetation_imaris.mat', 'imgseg')
   save('./Test/Images/Develop/Aryeh/140305_RUES2_36hBMP4_Bra_Snail_Sox2_surfaces_imaris_figure.mat', 'surfaces')
   
   %%
   save(savefile, 'img', 'imgseg', 'stats', 'surfaces')
end






%% Visualization / Plotting Statistics etc


%% Push Surfaces to imaris

if useimaris
   %%
   %nset = 10;
   nset = size(surf);
   sfset = surf(1:nset);
   fcset = fac(1:nset);
   nmset = norm(1:nset);

   imarissetsurface('Aryeh Segments', sfset, fcset, nmset);
end


%% Plotting Statistics as Colored Surfaces

if verbose
   %%
   figure(11)
   clf
   %implot3d(imgseg);
   %cdata = [stats.Centroid];
   %cdata = cdata(3,:);
   
   %cdata = [stats.MeanIntensity];
   cdata = [stats{3}.MeanIntensity];
   
   
   
   %cdata = [stats{2}.Volume];
   
  
   param.color.data = cdata;
   imsurfaceplot3d(imgseg, param)
end

%% Plotting Statistics as Scatter Plots
if verbose
   
   xmeasure = 'Volume';
   ymeasures = {'MeanIntensity', 'MinIntensity','MaxIntensity', 'MedianIntensity'};
   cccol = {'b', 'g', 'r', 'k'};
   
   for cc = 1:4
      figure(60+cc)
      clf
      set(gcf, 'Name', ['channel: ' num2str(cc)]);
      for i = 1:4
         ax(i) = subplot(2,2,i);
         plot([stats{cc}.(xmeasure)], [stats{cc}.(ymeasures{i})], ['*' cccol{cc}]);
         xlabel(xmeasure); ylabel(ymeasures{i});
      end
      linkaxes(ax, 'x')
   end  
end


%%
if verbose
   
   xmeasure = 'z';
   xdata = [stats{1}.Centroid];
   xdata = xdata(3,:);

   ymeasures = {'MeanIntensity', 'MinIntensity','MaxIntensity', 'MedianIntensity'};
   cccol = {'b', 'g', 'r', 'k'};

   
   for cc = 1:4
      figure(70+cc)
      clf
      set(gcf, 'Name', ['channel: ' num2str(cc)]);
      for i = 1:4
         ax(i) = subplot(2,2,i);
         plot(xdata, [stats{cc}.(ymeasures{i})], ['*' cccol{cc}]);
         xlabel(xmeasure); ylabel(ymeasures{i});
      end
      linkaxes(ax, 'x')
   end  
   
end


%end

