%% preprocessing

img = imread('./Test/Images/raw.tif');
imgd = double(img);
imgd = mat2gray(imgd);

imgdg = gaussianFilter(imgd,2,10);
imgmed = medianFilter(imgdg, 3);
imgmsf = meanShiftFilter(imgd, 3, 0.1);

imgdisk = diskFilter(imgmed, 3, 4, 1, -1);

imgmax = imregionalmax(imgdisk);


figure(50)
subplot(2,3,1)
hist(imgd(:), 256)
title('intensity')
subplot(2,3,2)
imglogvals = log2(imgd(:)+eps);
imglogvals(imglogvals < -15) = -15;
imglogvals(imglogvals > 0) = 0;
hist(imglogvals, 256)
title('log2 intensity')
subplot(2,3,3)
hist(log2(imgmed(:)+eps), 256)
title('log2 median filtered')
subplot(2,3,4)
plot(sort(imgd(:)))
subplot(2,3,5)
plot(sort(imglogvals))
subplot(2,3,6)
plot(sort(log2(imgmed(:)+eps)))

figure(51)
imshow([imoverlay(imgd, imgmax) imoverlay(imgdisk, imgmax)])



%% find maxima via log filter

%%% Laplacian of Gaussians
imglog = logFilter(1 - imgmed, 15);

%%% Detect Maxima
imgmax = imregionalmax(imglog);

%%% Combine nearby points and shrink to single points
%imgmax = imdilate(imgmax, strel('disk', 3));
%imgmax = bwmorph(imgmax,'shrink',inf);

%%% filter noise
maxvals = imgd(imgmax);
meanmaxvals = imfiltervalues(imgd, imgmax, 3);
threshold = mixtureOfGaussiansThreshold(log2(meanmaxvals+eps), 0.9)
threshold = 2^threshold

imgmax(imgmax) = (imgd(imgmax) >= threshold);
imgmaxlabel = bwlabel(imgmax);


%%% plotting the results
figure(51)
imshow([imoverlay(imgd, imgmax) imoverlay(imlog, imgmax)])

figure (52)
subplot(1,2,1)
hist(maxvals, 256)
subplot(1,2,2)
hist(log2(meanmaxvals+eps), 256)


%% ROI: Thresholding and Closing / Opening Original Image  

%thresholdEmpirical = 3000/ 65535.;

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


figure(54)
subplot(2,1,1)
imshow([imgd, imgthres])
subplot(2,1,2)
imshow([imoverlay(imgd, imgroi), imoverlay(imgthres, imgroi)])


%% segment using the maxima and watershed

% intensity method

% watershed on median filtered image

imgmin = imimposemin(1 - imgtmsf, imgmax);
imgws  = watershed(imgmin);

imgws = uint16(imgroi) .* imgws;

figure(53)
subplot(1,3,1)
imshow(imcolorize(imgws))
subplot(1,3,2)
imshow(double(imcolorize(imgws)) .* gray2rgb(imgd))
subplot(1,3,3)
imshow(imgmin)


%% segmentation using Propagation

distanceCutoff = Inf;

imgmedgrad = imgradient(imgmed);

mi = max(imgmed(:));
mg = max(imgmedgrad(:));

imgprop = imadd(imgmed, 5.0 * mi / mg * imgmedgrad);



[imgproplabel, dist] = segmentByPropagation(imgprop, imgmaxlabel, imgroi, 0.1, 0, 5);
imgproplabel(dist>distanceCutoff) = 0;

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


%% 

figure




%% segmentation using two step water shed





%% segmentation using water shed on gradient overlaid image

distanceCutoff = Inf;

imgmedgrad = imgradient(imgmed);

mi = max(imgmed(:));
mg = max(imgmedgrad(:));

imgprews = imsubtract(imgmed, .5 * mi / mg * imgmedgrad);

imgmin = imimposemin(1 - imgprews, imdilate(imgmax, strel('disk',3)));
imgws  = watershed(imgmin);

imgws = uint16(imgroi) .* imgws;

figure(53)
ax(1) = imsubplot(2,3,1);
imshow(imoverlay(imgmed, immax))
ax(2) = imsubplot(2,3,2);
imshow(imgmedgrad)
ax(3) = imsubplot(2,3,3);
imshow(imgprews);
ax(4) = imsubplot(2,3,4);
imshow(imgmin)
ax(5) = imsubplot(2,3,5);
imshow(imcolorize(imgws));
ax(6) = imsubplot(2,3,6);
imshow(imoverlay(double(imcolorize(imgws)) .* gray2rgb(imgd), immax));

%imshow(imerode(imgmedgrad, strel('disk', 1)))

linkaxes(ax, 'xy');

%% segmentation using ray method


figure(1)
imsubplot(1,2,1)
plot(linspace(1,10))

imshow(imgd)


%% stuff

imgmaxvals = df(imgmax);
threshold = mixtureOfGaussiansThreshold(imgmaxvals, 0.9);
%immax(immax > 0) = immaxvals >= threshold;


%montage([imoverlay(imgd, imgmax) imoverlay(imgmed, imgmax); imoverlay(imgmsf, imgmax) imoverlay(imgdisk, immax)])




%% other sutff


xy1 = [1, 0];
xy2 = [-1, -0.0001];

%xy1 = xy1 / norm(xy1);
%xy2 = xy2 / norm(xy2);

atan2(xy1(1) * xy2(2) - xy1(2) * xy2(1), dot(xy1, xy2)) + pi


%% testing

nrays = 15;
ds = 1;    % magnitude of change on ray
weight_angle = 1;
weight_grad = 1;


phi = 0:2*pi/nrays: (2*pi - 2*pi/nrays);
rays= [cos(phi); sin(phi)];

shape = ones(1,length(rays));
%shape(4) = 1;


% initialize
nrays = length(rays);
pos = [shape; shape] .* rays;

figure(32)
clf
impoly(gca, pos');


gradE = zeros(1,nrays);

diffvec = [pos(:, 2:end) pos(:, 1)] - pos;

dist2 = sum(diffvec.^2, 1);
dist = sqrt(dist2);
%normdiffvec = 1/dist;
%normdiffvec = [normdiffvec; normdiffvec] .* diffvec;


%%%
%%% engergy due to deviation form perfect polygon
%%%

% calculating the inner angles
dotv = sum([diffvec(:, end) diffvec(:, 1:end-1)] .* diffvec,1);
detv = [diffvec(1,end) diffvec(1,1:end-1)] .* diffvec(2,:) - [diffvec(2,end) diffvec(2,1:end-1)] .* diffvec(1,:);
angle = atan2(detv, dotv) + pi;  % outer angle
angle = 2 * pi - angle;                   % inner angle

% calculating the perturbed inner angles
pos_delta = pos + ds * rays;

diffvec_delta_lo = pos_delta - [pos(:, end) pos(:, 1:end-1)];
diffvec_delta_hi = [pos(:, 2:end) pos(:, 1)] - pos_delta;

dotv_delta = sum(diffvec_delta_lo .* diffvec_delta_hi,1);
detv_delta = diffvec_delta_lo(1,:) .* diffvec_delta_hi(2,:) - diffvec_delta_lo(2,:) .* diffvec_delta_hi(1,:);

angle_delta = arrayfun(@atan2, detv_delta, dotv_delta) + pi;  % outer angle
angle_delta = 2 * pi - angle_delta;                           % inner angle


% grad
angle_opt = (1-2/nrays) * pi % perfect n-gon inner angle

gradE = gradE + weight_angle * sin(angle - angle_opt).^2;


%%%
%%% surface tension
%%%

surface = sum(dist)
%surface_delta = 

%dist2_delta_lo = sum( diffvec_delta_lo.^2 ,1);
%dist2_delta_hi = sum( diffvec_delta_hi.^2, 1);
%dist2_lo = dist2;
%dist2_hi = [dist2(2:end-1)

%dist_delta = dist2_delta_lo + dist2_delta_hi + sum(dist2
% surface tension




% radial penalty


% area penalty


% intensity coverage


% gradient coverage



