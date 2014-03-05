%% Segmentation of 2D Images - Template
%
% Template to create your individual segmentation routine
%
% Save this file to ./Skripts/User/ under a descriptive name
% and modify to find the combination of steps that best segment your image
%
% Run this file cell by cell to experiment with segmentation parameters. 
% When satisfied comment out the display features and unwanted routines and
% create your segmentation function for processing all your images
%
% To help others, include in ./Test/Images a small sample image with 
% a unique and descriptive name.

%% Clean Start

clear all
close all
clc

initialize


%% Plotting

verbose = true;     % switch to turn on figure generation
figure_offset = 0;  % offset for figure numbering / use to compare results etc..


%% Load Image 
%
% result: - img    image of intensity values in [0,1]

filename = './Test/Images/StemCellNuclei.tif';

imgraw = imread(filename);
imgraw = double(imgraw);
imgraw = mat2gray(imgraw); % rescaling intensity range to lie within [0,1]

img = imgraw;

if verbose
   figure(1 + figure_offset)
   set(gcf, 'Name', ['Load: ' filename])
   imshow(img)
end

%% Preprocessing (optional)
%
% result: - img    image of intensity values in [0,1]

% appliy some filters to remove noise

% gaussian filter
%param.filter.ksize = [5, 5];    % size of the filter kernel (h x w box)
%param.filter.sigma = []         % std of gaussian [] = param.filter.ksize / 2 / sqrt(2 * log(2)); 
%img = gaussianFilter(img, param.filter.ksize, param.filter.sigma);

% mean shift filter - edge preserving 
%param.filter.ksize = 3;              % size of the filter kernel (h x w box)
%param.filter.intensity_width = 0.1;  % max deviaiton of intensity values to include in mean
%param.filter.iterations = 1;         % number of iterating the filtering 
%img = meanShiftFilter(img, param.filter.ksize, param.filter.intensity_width, param.filter.iterations);

% median filter - edge preserving
param.filter.ksize = 3;               % size of the filter kernel (h x w box) 
img = medianFilter(img, param.filter.ksize);

% bilateral filter - edge preserving
%param.filter.ksize = 3;              % size of the filter kernel (h x w box) 
%param.filter.sigma_space = [];       % std of gaussian in space  [] = param.filter.ksize / 2 / sqrt(2 * log(2));
%param.filter.sigma_intensity = [];   % std fo gaussian in intensity [] = 1.1 * std(img(:));
%img = bilateralFilter(img, param.filter.ksize, param.filter.sigma_space, param.filter.sigma_intensity);

% function filter
%param.filter.ksize = [5, 5];                    % size of the filter kernel (h x w box) 
%param.filter.function = @(x)(max(x,[],2)); ;    % function acting on array that for each pixel (1st dim) contains its neighbourhood in 2nd dim
%                                                % should return a vector of the new pixel values
%img = medianFilter(img, param.filter.ksize, param.filter.function);

% others: see ./Filtering folder

if verbose && total(abs(img-imgraw)) > eps
   figure(2 + figure_offset)
   set(gcf, 'Name', ['Preprocess: ' filename]);
   implottiling({imgraw, img},{'imgraw', 'img'});
end


%% Thresholding / Masking
%
% result: - imgmask    binary mask that specifies region of interest, get mask correct here
%         - imgth      thresholded image with removed background


% determine threshold using the histogram on logarithmic intensities

imgvalslog = log2(img(:)+eps);
imgvalslog(imgvalslog < -15) = -15; % bound lower values
imgvalslog(imgvalslog > 0) = 0;     % bound upper values

param.threshold.MoG = 0.5;          % probability of a pixel belonging to foreground
                                    % decrease to decrease the fitted threshold                        
%thlog = 2^thresholdMixtureOfGaussians(imgvalslog, param.threshold.MoG);  % this usually takes long


% direct thresholding methods (with optional prefiltering)

%param.filter.ksize = 3;
%imgf = medianFilter(img, param.filter.ksize);
imgf = img;

thotsu = thresholdOtsu(imgf);
thentropy = thresholdEntropy(imgf);
thmentropy = thresholdMutualEntropy(imgf);
thmog = thresholdMixtureOfGaussians(imgf, 0.5);


% thresholding using local maxima statistics

imgvalsmax = img(imregionalmax(img));
%imgvalsmax = img(imextendedmax(img, 0.01));
imgvalslogmax = log2(imgvalsmax);
param.threshold.MoG = 0.5;   
thmax = 2^thresholdMixtureOfGaussians(imgvalslogmax, param.threshold.MoG);


% select a threshold and create mask and thresholded image
th = thmax;

imgth = img;
imgth(imgth < th) = 0;
imgmask = imgth > 0;

% remove small fragments and small isolated minima with imopen
imgmask = imopen(imgmask, strel('disk', 2));   % larger disk size removes larger fragments
% close
%imgmask = imclose(imgmask, strel('disk', 2));  % larger disksize closes more
% dilate
%imgmask = imdilate(imgmask, strel('disk', 2)); % increase mask region 
% erode
%imgmask = imerode(imgmask, strel('disk', 2));  % decreasde mask region
% fill holes
%imgmask = imfill(imgmask,'holes');


if verbose

   prt = '\n\nthresholds:\n===========\n';
   thnames = {'thlog', 'thmog', 'thotsu', 'thentropy', 'thmentropy', 'thmax'};
   thdescription = {'MoG on log(img) = %7.5f', 'MoG(img)        = %7.5f', 'Otsu(img)       = %7.5f',... 
                    'Entropy on hist = %7.5f', 'Mutual entropy  = %7.5f', 'MoG on local max= %7.5f'};
   for t = 1:length(thnames)
      if exist(thnames{t}, 'var')
         prt = [prt '\n' sprintf(thdescription{t}, eval(thnames{t}))];
       end
   end
   prt = [prt '\n-------------------------\n'];
   prt = [prt sprintf('threshold       = %7.5f', th)];
   fprintf([prt '\n']);


   figure(10 + figure_offset)
   set(gcf, 'Name', ['Thresholding: ' filename])
   implottiling({imgth, imoverlay(img, imgmask)});
   
   figure(11 + figure_offset)
   set(gcf, 'Name', ['Thresholding Histograms: ' filename])   
   
   subplot(2,4,1);
   hist(imgraw(:), 256);
   title('raw intensities')
   subplot(2,4,5);
   plot(sort(imgraw(:))) % x-axis is effectively number of pixels <= ordinate in following plots
      
   if exist('imgvalslog', 'var')
      subplot(2,4,2);
      hist(imgvalslog, 256);
      title('log intensities');
      subplot(2,4,6);
      plot(sort(imgvalslog(:)))
   end
   if exist('imgvalsmax', 'var')
      subplot(2,4,3);
      hist(imgvalsmax, 256);
      title('local max intensities');
      subplot(2,4,7);
      plot(sort(imgvalsmax(:)))
   end
   if exist('imgvalslogmax', 'var')
      subplot(2,4,4);
      hist(imgvalslogmax, 256);
      title('local log max intensities');
      subplot(2,4,8);
      plot(sort(imgvalslogmax(:)))
   end
   
end   



%% Seeding
%
% result: - imgmax     binary image indicating seeds for segmentation


%%% optional pre filtering

imgf = imgraw;

% gaussian smoothing
%imgf = gaussianFilter(imgf,3,10);

% median filter / if note cumpted above or different parameter set
imgf = medianFilter(imgf, 3);

% mean shift 
%imgf = meanShiftFilter(imgf, 3, 0.1);


%%% center enhancing filter

% Laplacian of Gaussians (LoG) - more robust / use on inverse image !
%param.filter.ksize = [15, 15];       % size of the filter = diameter of nuclei
%param.filter.sigma = [];           % std of gaussian ([] = ksize / 4 / sqrt(2 * log(2)))
%imgf = logFilter(max(imgf(:)) - imgf, param.filter.ksize, param.filter.sigma);

% disk filter (consists of inner disk and optional outer ring to enhanve edges
%param.filter.ksize = 12;            % size of filer h or [h, w]
%param.filter.ring_width = 2;        % width fo ring (disk radius is determined by outer radius - ringh_width)
%param.filter.disk_weight = 1;       % weight on inner disk
%param.filter.ring_weight = -1;      % weight on outer ring
%imgf = diskFilter(imgf, param.filter.ksize, param.filter.ring_width, param.filter.disk_weight, param.filter.ring_weight);

% Difference of Gaussians (DoG) - similar to LoG but less robust
%param.filter.ksize = [15, 15];      % size of the filter
%param.filter.sigma_in = [];         % std of inner Gaussian ([] = 1/1.5 * sigma_out)
%param.filter.sigma_out = [];        % std of outer negative Gaussian ([] = ksize / 2 / sqrt(2 log(2)) )
%imgf = dogFilter(imgf, param.filter.ksize,  param.filter.sigma_in,  param.filter.sigma_out);

% sphere filter
param.filter.ksize = [5, 5];
ker = fspecial2('sphere', param.filter.ksize);
%ker = ker + fspecial2('disk',  param.filter.ksize, 1, 1, -1);
imgf = linearFilter(imgf, ker);


% normalize
imgf = mat2gray(imgf);


%%% Maxima detection

% h-max detection (only local maxima with height > hmax are considered as maxima
param.filter.hmax = 0.001;  %0.02;
imgmax = imextendedmax(mat2gray(imgf), param.filter.hmax);

%local max
%imgmax = imregionalmax(imgf);

% constrain to maxima within mask
imgmax = imgmax .* cast(imgmask, class(imgmax));

% Combine nearby points and shrink to single points
imgmax = imdilate(imgmax, strel('disk', 1));
% fill holes - combination of nearby points can lead to holes
imgmax = imfill(imgmax,'holes');

% imgmax = bwmorph(imgmax,'shrink',inf);           % extended regions usually give better segmentatrion results


% plot the results. Ideally have one seed per nucleus
if verbose  % CK I see a figure(20) and such further down is that intended
   
   figure(20 + figure_offset)
   set(gcf, 'Name', ['Seeding: ' filename])
   implottiling({imoverlay(imgraw, imgmax), imoverlay(imgf, imgmax)});

end

%% Postprocess Seeds

% Join seeds that can be connected by lines that does not cross boundary
% CK where is this used??

param.join.distance = 15;     % maximal distance of seeds to join 
param.join.threshold = -Inf;  % minimum intensity of seed to be considered for joining
param.join.distance = 15;     % maximal distance of seeds to join 
param.join.threshold = -Inf;  % minimum intensity of seed to be considered for joining



%% test stuff

% single point seeds

imgmaxp =  bwmorph(imgmax,'shrink',inf);
imglabel = bwlabeln(imgmaxp);

labdist = imlabeldistances(imgmaxp);


pos = immask2coords(imgmaxp);
pairs = imlabellocalpairs(imgmaxp, 15, 'dist');
np = length(pairs);  %CK ;



figure(43)
clf
ax(1) = imsubplot(1,2,1);
implot(imoverlay(img, imgmaxp))

hold on
for i = 1:np
   
   xy =[pos(pairs(i,1),:); pos(pairs(i,2),:)];
   line(xy(:,1), xy(:,2))
end

pairs(1,:);
pos(pairs(1,1),:);
pos(pairs(1,2),:);


%ax(2) = imsubplot(1,2,2);
%imshow(mat2gray(labdist))



% see if two points should be joined inlcude joining line into seed image






%% Segmentation on Image

% optional fitlering
%imgf = img;
imgf = medianFilter(img,3);


%dilating the maxima can improve segmentation
imgmaxd = imdilate(imgmax, strel('disk', 1));
%imgmaxd = imgmax;

% watershed
imgmin = imimposemin(max(imgf(:)) - imgf, imgmaxd);
imgws = watershed(imgmin);
imgseg = double(imgws) .* double(imgmask);

% watershed on gradient
imgmaxd = imdilate(imgmax, strel('disk', 3));
imggrad = imgradient(gaussianFilter(img, 3, 10) .* imgmask);
rm = imimposemin(imggrad, imgmaxd);
ws = watershed(rm);
imgsegg = double(ws) .* double(imgmask);


figure(20)
implottiling({imcolorize(imgseg), imcolorize(imgsegg)}, {'watershed on image', 'watershed image gradient'});


figure(21)
implottiling({imoverlaylabel(img, imgseg), double(imcolorize(imgseg)/255) .* double(gray2rgb(img))}, {[], 'watershed on img overlaid on img'})


%% Watershed Segmentation on Image + Gradient

%imgf = img;
imgf = medianFilter(img,3);

imgfgrad = imgradient(imgf);
mi = max(imgf(:));
mg = max(imgfgrad(:));
imgmix = imsubtract(imgf, 0.5 * (mi/mg) * imgfgrad);

imgmin = imimposemin(max(imgmix(:)) - imgmix, imdilate(imgmax, strel('disk',0)));
imgws = watershed(imgmin);
imgseg = double(imgws) .* double(imgmask);

figure(30)

implottiling({imoverlay(imgf, imgmax), imgfgrad, imgmix; ...
              imgmin, imcolorize(imgws), imoverlay(double(imcolorize(imgws)/255) .* gray2rgb(img), imgmax)})



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

% check the file for details
% edit ./Test/testRaySegmentation.m



%% Postprocess Segmentation and alternative diagnositcs

%clumb splitting


%


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
imgtmsf = meanShiftFilter(imgthres, 3, 0.1);


figure(100)
subplot(2,1,1)
imshow([imgd, imgthres])
subplot(2,1,2)
imshow([imoverlay(imgd, imgroi), imoverlay(imgthres, imgroi)])


ADD last statement below to the segmentation_2D.m

%% Clean up segmentation and alternative diagnositcs

 
lbl = bwlabel(imgseg > 0);
stats = regionprops(logical(lbl), 'Area', 'PixelIdxList');
min_area = 40;
keep = [stats.Area] >= min_area;
% remove small nucs from lbl
for i = find(~keep)
    lbl(stats(i).PixelIdxList) = 0;
end
lbl = imfill(lbl, 'holes');



%% Extended Maxima to Join Label

imgf = mat2gray(medianFilter(imgraw,3));

imgfgrad = imgradient(imgf);
mi = max(imgf(:));
mg = max(imgfgrad(:));
%imgf = imsubtract(imgf, 0.5 * (mi/mg) * imgfgrad);


imglabel = bwlabeln(imgmax);
label = imlabel(imglabel);

% in each labeled region find max and set all values to this max
for l = label
   idx = find(imglabel == l);
   vals = imgf(idx);
   imgf(idx) = max(vals);
end

% use log or disk fitler again ??

imgff = imgf;
imgff = medianFilter(imgf,3);
%imgff = logFilter(max(imgf(:)) - imgf, [15,15]);
%imgff = diskFilter(imgf, 15, 1, 1, -1);
imgff = mat2gray(imgff);


% now find extended maxima

imgmm = imgff;
%imgmm = imhmin(imgmm, 0.05);
imgmm = imextendedmax(imgmm, 0.01);
imgmm = imdilate(imgmm, strel('disk', 0));
imgmm = imgmm .* cast(imgmask, class(imgmm));




figure(50)
clf
ax(1) = imsubplot(1,3,1);
imshow(imgff)
ax(2) = imsubplot(1,3,2);
imshow(imoverlaylabel(imgf, imglabel))
ax(3) = imsubplot(1,3,3);
imshow(imoverlay(imgraw, imgmm))

linkaxes(ax, 'xy');



%%
clf
ker = fspecial3('disk', 3, 0, 1, -1);

ker = fspecial2('dog', 15)

ker = fspecial2('sphere', 15)
%implot3d(mat2gray(ker))
%immontage(mat2gray(ker))
imshow(mat2gray(ker))

size(ker)

%%
figure

imshow(imgmask)



