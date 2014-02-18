%% Testing 3D Segmentation

% conventions
% figure numbering:
%  1-10 raw data loading
% 11-20 filtering / preprocessing
% 21-30 maxima detection
% 31-40 segmenation
%


%% clear 

close all
clear all
clc

fp = [1   705   560   420];
set(0, 'DefaultFigurePosition', fp)


%% synthetic data
[stack, ctr] = syntheticConfocalData();

figure(1)
clf
imshow3d(stack)

figure(2)
immontage(stack)


%% test data - "noisy"
img= imread('./Test/Images/130905_mES_TGFLEFGFP/130905_mES_TCFLEFGFP_t001_z006_c001.tif');
figure
image(img)

clear all
close all
clc

x0 = 170; y0 = 320; wh = 70;
stack = imreadstack('./Test/Images/130905_mES_TGFLEFGFP/*_t001*c002.tif', 'PixelRegion', {[x0 x0+wh], [y0 y0+wh]});

figure(1)
clf
imshow3d(stack, 'BoxRatios', [1,1,0.2] )
xlabel('x'); ylabel('y'); zlabel('z');

figure(2)
immontage(max(stack(:)) - stack)


%%

stack = syntheticConfocalData();



%% thresholding 

figure(3)
hist(log2(stack(stack > 0)), 120)

%th= 2^thresholdMixtureOfGaussians(log2(stack(stack>0)), 0.1)

th = 1;

figure(4)
clf
stackth = stack;
stackth(stack < th) = 0;
imshow3d(stackth)


%% clear 
clear all
close all
clc

%% other test data - "intensity"

x0 = 180; y0 = 210; wh = 100;
stack = imreadstack('./Test/Images/20130911T114814/W1F127T0021Z*C1.tif', 'PixelRegion', {[x0 x0+wh], [y0 y0+wh]});

figure(1)
clf
imshow3d(stack, 'BoxRatios', [1,1,0.2] )
xlabel('x'); ylabel('y'); zlabel('z');

figure(2)
immontage(stack)

figure(3)
hist(log2(stack(stack > 0)), 120)

th = 2^thresholdMixtureOfGaussians(log2(stack(stack>0)), 0.1);


figure(4)
clf
stackth = stack;
stackth(stack < th) = 0;
imshow3d(stackth)

%%

figure(5)
clf
stackmed = medianFilter(stackth, [3 3 1]);
imshow3d(stackmed)

figure(6)
immontage(stackmed)




%% detect 3d log maxima 

stacklog = logFilter(max(stack(:)) - stack, [20 20 17]);
stackmax = imregionalmax(stacklog);

figure(104)
clf
colormap(autumn)
is = stacklog- mean(stacklog(:));
is(is < 0) = 0;
imshow3d(is, 'BoxRatios', [1,1,0.2]);


figure(95)
clf
imshow3d(stackmax , 'BoxRatios', [1,1,0.2] )

% noise

ker = ones(3, 3, 3); ker = ker / sum(ker(:));
iv = imfiltervalues(stack, stackmax, ker);

%th = thresholdMixtureOfGaussians(iv, 0.)
th = thresholdOtsu(iv)
th = thresholdMutualEntropy(iv)
th = thresholdEntropy(iv)
th = 1;

pos = find(stackmax);

stackmaxf = stackmax;
stackmaxf(pos(iv < th)) = 0;

figure(24)
hist(iv, 150);

figure(25)
clf
imshow3d(stackmaxf)

sum(stackmaxf(:))

%% 3d water shed 

mask = stack>0;
rm = imimposemin(max(stack(:)) - stack, stackmaxf);
ws = watershed(rm);


figure(97)
clf
imshow3d(double(ws) .* double(mask))

%length(unique(ws .* uint16(mask)))-1


%% check active contours

img = mat2gray(stack(:,:,8));

figure(10)
imshow(img)


imgd = double(mat2gray(img));

imgdg = gaussianFilter(imgd,2,10);
imgmed = medianFilter(imgdg, 3);
imgmsf = meanShiftFilter(imgdg, 3, 0.05);

imglog = logFilter(1 - imgmed, 8);

%%% Detect Maxima
imgmax = imregionalmax(imglog);

iv = imfiltervalues(imgmed, imgmax, 3);
th = 2^thresholdMixtureOfGaussians(log2(iv));

mask = imgmed > th;

pos = find(imgmax);
imgmax(pos(iv<=th)) = 0;


figure(11)
clf
ax(1) = imsubplot(3,3,1);
imshow(imoverlay(imgd, imgmax))
ax(2) = imsubplot(3,3,2);
imshow(imgdg)
ax(3) = imsubplot(3,3,3);
imshow(imgmed)
ax(4) = imsubplot(3,3,4);
imshow(imgmsf)
ax(5) = imsubplot(3,3,5);
imshow(imoverlay(imglog, imgmax))

imsubplot(3,3,6);
hist(log2(iv), 128)

imsubplot(3,3,7)
imshow(mask)

linkaxes(ax, 'xy')






%% test segmentation

pos = find(imgmax, 1, 'first');
start = zeros(h,w);
start(pos) = 1;
start = imdilate(start, strel('disk', 2));

%seg = activecontour(imgmed, start, 'Chan-Vese'); %'Chan-Vese' 'edge'

ws = imimposemin(imcomplement(img), imgmax);
ws = double(watershed(ws)) .* mask;

figure(12)
ax(1) = imsubplot(2,2,1);
imshow(img)
ax(2) = imsubplot(2,2,2);
imshow(imoverlay(imgmsf, imgmax));
ax(3) = imsubplot(2,2,3);
%imshow(imoverlay(imgmed, seg));
imshow(double(imcolorize(ws)) .* gray2rgb(imgd))
ax(4) = imsubplot(2,2,4);
imshow(imcolorize(ws));

linkaxes(ax, 'xy')





%% segement 2d images

[h,w,l] = size(stack);
seg = zeros(h,w,l);
imax = zeros(h,w,l);

show = true;

for zi = 1:l
%for zi = 6:6
   img = stack(:,:,zi);
   imgd = double(mat2gray(img));

   imgdg = gaussianFilter(imgd,2,10);
   imgmed = medianFilter(imgdg, 3);
   imgmsf = meanShiftFilter(imgd, 3, 0.1);

   %imgdisk = diskFilter(imgmed, 3, 4, 1, -1);
   %imgmax = imregionalmax(imgdisk);

   %%% Find Maxima vial log filter 
   %%% Laplacian of Gaussians
   imglog = logFilter(1 - imgmed, 4);

   %%% Detect Maxima
   imgmax = imregionalmax(imglog);

   %%% Combine nearby points and shrink to single points
   %imgmax = imdilate(imgmax, strel('disk', 4));
   %imgmax = bwmorph(imgmax,'shrink',inf);

   %%% filter noise
   maxvals = imgd(imgmax);
   meanmaxvals = imfiltervalues(imgd, imgmax, 3);
   threshold = thresholdMixtureOfGaussians(log2(meanmaxvals+eps), 0.9);
   threshold = 2^threshold

   threshold = 0.01;
   
   imgmax(imgmax) = (imgd(imgmax) >= threshold);
   imgmaxlabel = bwlabel(imgmax);


   if show
   figure(51)
   imshow([imoverlay(imgd, imgmax) imoverlay(imglog, imgmax)])

   figure (52)
   subplot(1,2,1)
   hist(maxvals, 256)
   subplot(1,2,2)
   hist(log2(meanmaxvals+eps), 256)
   end

   %thresholdEmpirical = 3000/ 65535.;

   imgroi = im2bw(imgmed, threshold);


   % larger sizes exclude dividing cells !!
   imgroi = imerode(imgroi, strel('disk',2));
   imgroi = imdilate(imgroi, strel('disk',2));
   %imgroi = bwmorph(imgroi, 'open', 15);
   %imgroi = bwmorph(imgroi, 'close');

   imgthres = zeros(size(imgd));
   imgthres(imgroi) = imgd(imgroi);

   imgtmed = medianFilter(imgthres,3);
   imgtmsf = meanShiftFilter(imgthres, 3, 0.1);

   if show
   figure(54)
   subplot(2,1,1)
   imshow([imgd, imgthres])
   subplot(2,1,2)
   imshow([imoverlay(imgd, imgroi), imoverlay(imgthres, imgroi)])
   end
   
   
   %%% segmentation using Propagation
   distanceCutoff = Inf;

   imgmedgrad = imgradient(imgmed);

   mi = max(imgmed(:));
   mg = max(imgmedgrad(:));

   imgprop = imadd(imgmed, 5.0 * mi / mg * imgmedgrad);


   [imgproplabel, dist] = segmentByPropagation(imgprop, imgmaxlabel, imgroi, 0.1, 0, 5);
   imgproplabel(dist>distanceCutoff) = 0;

   if show
   figure(53)
   ax(1) = subplot(1,3,1);
   imshow(imoverlay(imcolorize(imgproplabel), imgmax, [1 1 1]))
   ax(2) = subplot(1,3,2);
   imshow(double(imcolorize(imgproplabel)) .* gray2rgb(imgd))

   distinf0 = dist;
   distinf0(dist==Inf) = 0;
   ax(3) = subplot(1,3,3);
   imshow(mat2gray(distinf0))

   linkaxes(ax, 'xy')

   figure(54)
   imshow([imgmedgrad, imgprop])
   end


   % save the segments
   seg(:,:,zi) = imgproplabel;
   imax(:,:,zi) = imgmax;
   
end


figure(60)
imshow3d(seg)


figure(61)
ov = zeros(h,w,3,l);
for z = 1:l
   ov(:,:,:,z) = imoverlay(seg(:,:,zi), imax(:,:,zi));
end
montage(ov)

%% identify 3D segments form 2D

segE = zeros(h,w,l);
for zi=1:size(seg,3)
   segE(:,:,zi) = imerode(seg(:,:,zi)>0, strel('disk',2)) .* seg(:,:,zi);
end

label = bwlabeln(seg > 0.5);

figure(61)
clf
imshow3d(label)



%% combine 2d segments in previous component

nlabel = max(label(:))

mask = (label == 3);
seglabel = seg .* mask;

%[minp, maxp] = imboundingbox(seglabel);
[seglabel, minp, maxp] = imextract(seglabel, 'BoundingBox');

figure(62)
clf
imshow3d(seglabel);


%% 


% parameter
min_overlap = 5; 


nz = size(label,3);
labs = unique(seglabel(:));
labs = labs(labs~=0)

cells = [];

slice = seglabel(:,:,1);

% in first slice all segments will be labeled as seeds for cells
cells = unique(slice);
cells = cells(cells~=0);
 
for zi = 2:nz
   
   slice_next = seglabel(:,:,zi);
   label_next = unique(slice_next);
   label_next = label_next(label_next~=0);
   
   % compute overlaps of seeds with labels in next slice

   nl = length(label_next);
   for l = nl:-1:1
      mask_next(l, :, :) = (slice_next == labels_next(l));
      area_next(l) = sum(sum( mask_next(l, :, :) ));
   end
   nc = length(cells);
   for c = nc:-1:1
      mask_cells(c, :, :) = (slice == cell(c));
      area_cells(c) =  sum(sum( mask_cells(c, :, :) ));
   end   
   
   clear overlap
   for c = nc:-1:1
      for l = nl:-1:1
         overlap(c, l) = sum(sum(mask_cells(c, :, :) .* mask_next(l, :,:)));
      end
   end
   
   
   for l=1:nl
      assigned  = 0;
      
      for i = ind
         
         ovlp = overlap(c,l);
         
         % next segment covers the cell entirely
         % could check for to large coverage -> indicating splitting could use this for optimal tracking
         if ovlp == area_cells(c) && ovlp <= area_next(l) 
            slice_next(mask_next(l)) = c; % assinge cell label to next segment 
            assigned = assigned + 1;
         elseif ovlp             
            
            
         end
      end
      
      if assinged == 0 % could not assign label -> new seed
         cells = [cells label_next(l)];
      end
      
   end
   
   
   
      
   
   
   




% move through z stack and decide on new merging new cells etc




