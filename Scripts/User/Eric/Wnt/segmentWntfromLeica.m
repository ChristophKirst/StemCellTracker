function imgseg = segmentWntfromLeica(img)

imgd = double(img);
imgd = mat2gray(imgd);

imglogvals = log2(imgd(:)+eps);
% ie lower limit 1 on 16 bit image.
imglogvals(imglogvals < -15) = -15;
imglogvals(imglogvals > 0) = 0;

%% thresholding / masking 

param.filter.median = 5; %3;
imgmed = medianFilter(imgd, param.filter.median);

thlog = 2^thresholdMixtureOfGaussians(imglogvals, 0.5);
thmed = thresholdMixtureOfGaussians(imgmed, 0.5);

imgth = imgd;
imgth(imgth < thmed) = 0;

imgmask = imgth > 0;
% remove small fragments and small isolated minima with imopen
imgmask = imopen(imgmask, strel('disk', 2));


figure(1)
implottiling({imgth, imoverlay(imgth, imgmask)})


imgf = immask(imgmed, imgmask);
imgf = medianFilter(imgf, 3);
imgf = immask(imgf, imgmask);
figure(2)
implot(imgf)

%% Filtering / Seeding

% gaussian smoothing
%imgdg = gaussianFilter(imgd,3,10);
%imgdg = img;
imgdg = mat2gray(imclip(mat2gray(img), 0, 0.4));

figure(7)
imshow(imgdg)

% median filter  ?? also computed above
% imgf = medianFilter(imgdg, 3);

% mean shift 
%imgf = meanShiftFilter(imgd, 3, 0.1);

% disk [ksize or outer-box, width of ring, wt-inner disk, wt-ring]
%imgf = diskFilter(imgf, 11, 1, 1, -1);
%imgf = sphereFilter(imgdg, 11);

% Laplacian of Gaussians
param.filter.logsize = [20, 20];
imgf = logFilter(max(imgdg(:)) - imgdg, param.filter.logsize);


% find maxima using h-max detection 
param.filter.hmax = 0.02;  %0.02;
%imgmax = imregionalmax(imgf);
imgmax = imextendedmax(mat2gray(imgf), param.filter.hmax);


figure(6)

implottiling({imoverlaylabel(mat2gray(imgf), bwlabeln(imgmax)), imoverlaylabel(mat2gray(img), bwlabeln(imgmax))});


%% Watershed Segmentation on Image

%dilating the maxima can help
imgmaxd = imdilate(imgmax, strel('disk', 1));
%imgmaxd = imgmax;

% watershed
imgmin = imimposemin(max(imgmed(:)) - imgmed, imgmaxd);
imgws = watershed(imgmin);
imgseg = double(imgws) .* double(imgmask);

imgseg = imlabelseparate(imgseg);

% %% Clean up segmentation and alternative diagnositcs
% 
% lbl = bwlabel(imgseg > 0);
% stats = regionprops(logical(lbl), 'Area', 'PixelIdxList');
% min_area = 40;
% keep = [stats.Area] >= min_area;
% % remove small nucs from lbl
% for i = find(~keep)
%     lbl(stats(i).PixelIdxList) = 0;
% end
% lbl = imfill(lbl, 'holes');
% 
% [bndry, lbl] = bwboundaries(lbl > 0, 8, 'noholes');
% stats = regionprops(logical(lbl), 'Area', 'Centroid', 'PixelIdxList');
% imgmask = lbl > 0;

%%

figure(11)
clf
implottiling({imcolorize(imgseg), gray2rgb(mat2gray(img)) .* imcolorize(imgseg)})
