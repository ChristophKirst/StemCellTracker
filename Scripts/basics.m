%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Basic Setup and Processing Steps %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

initialize

%% Generate Experiment


exp = Experiment('Name', 'Example', 'Date', datestr(now), 'Description', 'Test the Code',...
                 'BaseDirectory', './Test/Data/Experiment', ...
                 'ImageDirectoryName', '../../Images/hESCells_tif_stack',...
                 'ResultDirectoryName', '.', ...
                 'ReadImageCommandFormat', 'imread(''<directory>/<file>'')',...
                 'ImageFileFormat', 'W1F127T<time,4>Z<z,2>C1.tif', ...
                 'ImageFormatNames', {'z', 'time'},...
                 'ImageFormatRange', {1:15, 1});

exp.Info()

%% Save/Load some test data
 
sdata = rand(1,3);

exp.SaveResult('testdata.mat', sdata)

ldata = exp.LoadResult('testdata.mat');
ldata - sdata

%% Load Images an plot

exp.ReadImageCommand([6, 1])

img = exp.ReadImage([6, 1]);
size(img)

figure(1); clf;
implot(img)

%% Read an Image Stack and plot in 3d

for z = exp.ImageFormatRange{1}
   img(:,:, z) = exp.ReadImage([z, 1]);
end
size(img)

imgr = imresample(img, [3, 3, 2]);
size(imgr)

figure(2); clf
implot3d(mat2gray(imgr));


%% Simple 2d Segmentation

% see other example scritps in ./Scrtips folder for more details, e.g.
% for 2d segmentation click on segment2D.m and press Strg + d
% for 3d segmentation click on segment3D.m and press Strg + d
% for tracking click on tracking.m and press Strg + d 

%loading
img = exp.ReadImage([6, 1]);
size(img)

%% Masking
imgf = mat2gray(log(double(img)));
%imgf = medianFilter(imgf, 2);
imgf = meanShiftFilter(imgf, 3);
imgf = mat2gray(imgf);
imgf = imclip(imgf, 0.1, 0.85);
imgf = mat2gray(imgf);

%th = thresholdMixtureOfGaussians(imgf)
th = 0.2;
imgmask = imgf >= th;
imgmask = imopen(imgmask, strel('disk', 3));

imgth = immask(imgf, imgmask);

imggrd = mat2gray(imgradient(imgth));
%imggrd = imclose(imggrd, strel('disk', 2));
imggrd = imclip(imggrd, 0, 0.25);
imggrd = mat2gray(imggrd);

figure(1)
implottiling({img, imgf; imoverlay(imgf, imgmask), imggrd})


%% Seeding

%imgs = diskFilter(imgth, 20, 1, 1, -2);
imgs = mat2gray(imgth) - 0.5 * imggrd;
imgs(imgs<0) = 0;
imgs = mat2gray(imgs);
imgs = imopen(imgs, strel('disk', 2));

%imgs = medianFilter(imgs, 2);

%imgs2 = sphereFilter(imgs, 10);
%imgs2 = diskFilter(imgs, 5);
imgs2 = logFilter(max(imgs(:))- imgs, 20);
imgs2 = mat2gray(imgs2);

imgmax = imextendedmax(imgs2, 0.05);
imgmax = imdilate(imgmax, strel('disk', 2));
imgmax = immask(imgmax, imgmask);

figure(13);
implottiling({imggrd, imgs, imgs2;  imoverlay(img, imgmax), imoverlay(img, imgmax), imoverlay(imgf, imgmax)})

%% Watershed

imgw = imgf;

imgmin = imimposemin(max(imgw(:)) - imgw, imgmax);
imgws = watershed(imgmin);
imgseg = immask(imgws, imgmask);

figure(30)
colormap jet
implottiling({imoverlay(imgw, imgmax), imcolorize(imgseg), imoverlaylabel(img, imgseg, true)})

%% Postprocess and create intial statistics

param = setParameter('volume.min',    50,...     % minimal volume to keep (0)
                     'volume.max',    inf,...    % maximal volume to keep (Inf)
                     'intensity.min', -inf, ...  % minimal mean intensity to keep (-Inf)
                     'intensity.max', inf, ...   % maximal mean intensity to keep (Inf)
                     'boundaries',    false, ... % clear objects on x,y boundaries (false)
                     'fillholes',     true,...   % fill holes in each z slice after processing segments (true)
                     'relabel',       true);    % relabel from 1:nlabelnew (true)

[imgpost, stats] = postProcessSegments(imgseg, param);

stats

figure(31)
colormap jet
implottiling({imoverlay(img, imgmax), imcolorize(imgpost), imoverlaylabel(img, imgpost, true)})


%% Create Objects and Frame form labeled image

param = setParameter('time' ,  0, ...   % time for objects (0)
                     'rescale',1, ...   % rescale coordinates r by this factor ([1, 1(, 1)])
                     'method', 'median'); % how to calcualte the intensity in Object, a string of any function, 'none' = dont calcualte ('median')

objs = label2Objects(imgpost, img, stats, param);

frame = Frame('objects', objs, 't', 0);

exp.Result = frame;


%% plot some statistics

figure(42)
hist(double([frame(1).objects.intensity]), 255)

%%

exp.Info()





