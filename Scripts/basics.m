%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Besic Setup of an Experiment and Functionality %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%% Load Images

exp.ReadImageCommand([6, 1])

img = exp.ReadImage([6, 1]);
size(img)

figure(1); clf;
implot(img)


%% Read the full Stack and plot

for z = exp.ImageFormatRange{1}
   img(:,:, z) = exp.ReadImage([z, 1]);
end
size(img)

imgr = imresample(img, [3, 3, 2]);
size(imgr)

figure(2); clf
implot3d(mat2gray(imgr));


%% Simple 3D Segmentation

% see other example scritps in ./Scrtips folder for more details and full functionality

%masking
img = exp.ReadImage([6, 1]);
size(img)

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

%% Statistics

stats = 






%%









