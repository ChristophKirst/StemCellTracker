%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Segment Diluted Cells in 2d %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clear classes
close all
clc

initialize

%% File Handler
fh = FileHandler('BaseDirectory',          './Test/Data/Experiment', ...
                 'ImageDirectoryName',     '../../Images/hESCells_Cytoo',...
                 'ResultDirectoryName',    '.', ...
                 'ReadImageCommandFormat', 'imread(''<file>'')',...
                 'ReadImageFileFormat',    'Big1_CFP_<time>.tif');

fh.info()

%%
fh.ReadImageCommand('time', 1)
fh.fileName('time', 1)
fh.fileName

%%
img = fh.readImage('time', 3);
figure(42); clf; colormap gray
implot(img)


%%
img = fh.readImage('time', 3);
img = mat2gray(img);

figure(1)
colormap(gray)
implot(img)

figure(2)
subplot(1,2,1);
hist(img(:), 256)
subplot(1,2,2);
hist(log(img(:))+eps, 256)

%%

th = 2^-2.8
imgth = img;
imgth(imgth < th) = th;
imgth = mat2gray(imgth);

figure(3)
colormap gray
implot(imgth)


imgmask = img > th;
%imgmask = imclose(imgmask, strel('disk', 1));
imgmask = imopen(imgmask, strel('disk', 5));

imglab = bwlabeln(imgmask);
imglab = postProcessSegments(imglab, setParameter('volume.min', 150, 'relabel' , false));
imgmask = imglab > 0;

figure(4)
colormap gray
implottiling({img, imoverlay(imgth, imgmask, 'r',  false)})



%% Circular detection using Houh transform
[centers, radii, metric] = imfindcircles(imgth,[3 15]);

figure(12)
imshow(imgth)
viscircles(centers, radii, 'EdgeColor','b');

centers = floor(centers(:,[2,1]));
%centers(centers < 1) = 1;


imgc = zeros(size(img));
idx = sub2ind(size(img), centers(:,1), centers(:,2));
imgc(idx) = 1;
imgc = imdilate(imgc, strel('disk', 3));

%
figure(5); clf; colormap jet
implot(imoverlaylabel(img, imgc))


%% Line Joining

imglab = bwlabeln(imgc);

param = setParameter('threshold.min',        0.1, ... % if profile comes below this absolute intensity objects are different
                     'threshold.max',        0.9,...    % if profile is always above this threshold, objects are joined
                     'threshold.change'    , 0.3, ...   % maximal rel change in intensitiy above objects are assumed to be different
                     'threshold.gradient'  , 0.4, ...   % maximal absolute gradient change below objects are joined, only if gradient image is supplied
                     'cutoff.distance'     , 30, ...    % maximal distance between labels (= 20)
                     'averaging.ksize'     , 3, ...     % ksize to calculate reference mean intensity (=3)
                     'addline'             , true);     % add a line between joined label (true)

[imgjoin, pairs, joins] = joinSeedsByRays(imglab, imgenh, param);


figure(3); clf; colormap jet
implot(imoverlaylabel(imgenh, imgjoin))

figure(4); clf; colormap jet
imsubplot(1,2,1); 
implot(imoverlay(img, imglab));
for i = 1:size(pairs, 1);
   col = hsv2rgb([i/size(pairs,1), 1, 1]);
   plotSeedPairs(imglab, pairs(i, :), col)
end

imsubplot(1,2,2)
implot(imoverlay(img, imglab));
for i = 1:size(joins,1)
   %col = hsv2rgb([i/size(joins,1), 1, 1])
   plotSeedPairs(imglab, joins(i,:),  col)
end




%% Doesnt work -> try simple stuff ignore cluster

imgs = imgth;
imgs = medianFilter(imgs);

imgs = diskFilter(imgs, [15, 15], 1, 1, 0);
imgs(imgs< 0) = 0;
imgs = mat2gray(imgs);

imgmax = imextendedmax(imgs, 0.01);


figure(5)
colormap gray
implottiling({imoverlay(imgth, imgmax, 'r', true), imoverlay(imgs, imgmax, 'r', true)})


%% Watershed

imgw = imgth;

imgmin = imimposemin(max(imgw(:)) - imgw, imgmax);
imgws = watershed(imgmin);
imgseg = immask(imgws, imgmask);
imgseg = imrelabel(imgseg);

figure(30)
colormap jet
implottiling({imoverlay(imgw, imgmax), imcolorize(imgseg), imoverlaylabel(img, imgseg, true)});

%% Postprocess and create intial statistics

param = setParameter('volume.min',    0,...     % minimal volume to keep (0)
                     'volume.max',    inf,...    % maximal volume to keep (Inf)
                     'intensity.min', -inf, ...  % minimal mean intensity to keep (-Inf)
                     'intensity.max', inf, ...   % maximal mean intensity to keep (Inf)
                     'boundaries',    false, ... % clear objects on x,y boundaries (false)
                     'fillholes',     true,...   % fill holes in each z slice after processing segments (true)
                     'relabel',       true);    % relabel from 1:nlabelnew (true)

[imgpost, stats] = postProcessSegments(imgseg, param);


figure(31)
colormap jet
implottiling({imoverlay(img, imgmax), imcolorize(imgpost), imoverlaylabel(img, imgpost, true)});

max(imgseg(:))
max(imgpost(:))


%% Create Objects and Frame form labeled image

param = setParameter('time' ,  0, ...   % time for objects (0)
                     'rescale',1, ...   % rescale coordinates r by this factor ([1, 1(, 1)])
                     'method', 'median'); % how to calcualte the intensity in Object, a string of any function, 'none' = dont calcualte ('median')

objs = label2Objects(imgpost, img, stats, param);

frame = Frame('objects', objs, 't', 0);

exp.Result = frame;


%% plot some statistics

figure(42); clf; 
subplot(1,2,1)
hist(double([frame(1).objects.intensity]), 15)
subplot(1,2,2)
hist(double([frame(1).objects.volume]), 15)

















%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Other Tries potentially usefull %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Enhance Image by Gradient
imggrad = imgradient(img);
imggrad = mat2gray(imtophat(imggrad, strel('disk', 10)));
%imggrad = mat2gray(imdilate(imggrad, strel('disk', 1)));

imggrad(imggrad > 0.5) = 0.5;
imggrad = mat2gray(imggrad);
imggrad = immask(imggrad, imgmask);
%imggrad = imopen(imggrad, strel('disk', 1));

%imgedg = edge(imggrad, 'log', 0.005);
imgenh = mat2gray(img) - imggrad;
imgenh(imgenh < 0.1) = 0.1;
imgenh = mat2gray(imgenh);


figure(5)
hist(imgenh(:),256)

%
%imgenh = mat2gray(img)  -0.5* mat2gray(imggrad);
%imgenh(imgenh< 0.05) = 0;
%imgenh = immask(imgenh, imgmask);
%imgenh = functionFilter(imgenh, [3,3], 'median');


figure(4);
implottiling({img, imggrad, imgenh})


%%

ws = watershed(imhmin(iminvert(imggrad), 0.5));
ws = immask(ws, imgmask);

figure(7)
implot(imcolorize(ws))


%%

%imgs = functionFilter(immask(imgenh, imgmask), [3,3], 'min');
imgs = diskFilter(imgenh,  [7,7], 0, 2, -0.2);
imgs = mat2gray(imgs);

imgmax = imextendedmax(imgs, 0.1);


figure(9); clf; colormap gray
implottiling({imoverlay(img, imgmax), imgs})






%% Distance transform 

imgdist = bwdist(iminvert(imgmask));
imgmax = imextendedmax(imgdist, 4);

% skeletonization
imgske = bwmorph(imger, 'skel', inf);
for i = 1:6
   imgend = bwmorph(imgske, 'endpoints');
   imgske = imgske - imgend;
   imgske(imgske < 0) = 0;
end

figure(5); clf
implottiling({imgdist, imgmax, imgske, diskFilter(imgske, [5,5])})

%% Ring Dilation

rout = 6;
rin = 5;

se = strel('disk', rout);
se = se.getnhood;
se2 = strel('disk', rin);
se2 = padarray(se2.getnhood, (rout-rin) * [1,1],0);
ser = se-se2;
ser = strel('arbitrary', ser);

imger = imdilate(imgmask, ser);
imger = imerode(imgmask, ser);
figure(15)
implot(imoverlay(img, imger))


%%

imggrad = imgradient(medianFilter(imgth, 3));
imgrad = mat2gray(imggrad);
%imggrad = edge(imgth, 'log', 0, 1);

%imggrad = imclose(imggrad, strel('disk', 1));


figure(13)
colormap gray
implottiling({imgth, imggrad})

%%

imgs = imgth - 0.0 * imggrad;
imgs(imgs<0) = 0;

figure(14)
colormap gray
implot(imgs)



