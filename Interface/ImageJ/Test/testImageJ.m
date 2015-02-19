%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test ImageJ interface %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initialize

%% Clean Classes

clear all
close all
javarmdynamicclasspath

%% Compile

ijcompile(1)


%% ImageJ Path

ijpath

%% ijinitialize

ijinitialize

%% Test 3D Image

img = syntheticLabeledImage([30, 50, 10], 20, 10);
figure(1); clf
colormap gray
implot3d(img)
freezecolormap

%% ijplot3d
ijplot3d(img, 'PixelDepth', 2)

%% ijplot3d color

figure(2)
imgc = imcolorize(img);
implot3d(imgc)

ijplot3d(imgc, 'PixelDepth', 2)

%% ijplot5d

imgmovie = repmat(imgc, 1, 1, 1, 1, 10);
size(imgmovie)

ijplot5d(imgmovie)


%% ijplot5d large data

img = syntheticLabeledImage([512, 512, 20], [50, 50, 4], 60);
imgc = imcolorize(img);

%min(imgc(:))
%max(imgc(:))

tic
imgmovie = repmat(imgc, 1, 1, 1, 1, 10);
ijplot5d(imgmovie, 'class', 'single')
toc

%%
tic
imgmovie = repmat(imgc, 1, 1, 1, 1, 10);
ijplot5d(imgmovie, 'class', 'uint8')
toc

