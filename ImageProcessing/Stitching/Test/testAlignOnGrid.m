%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test Image Alignment on Grid %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
clear classes
close all

initialize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Images On a Pre-aligned Grid

% note: convetion for girds is grid{x,y,z}
%       consistent with pql coordinates usage


%% Generate 2x2 Grid

img = loadTestImage();
img = mat2gray(img);

size(img)

figure(1);
implot(img);

imgs{1,1} = img(1:230,1:250);
imgs{2,1} = img(200:end,20:260);
imgs{1,2} = img(10:240,230:end);
imgs{2,2} = img(230:end,240:end);

figure(2);
implottiling(imgs, 'link', false);


%% Prepare image alignment

[img1, img2, per] = align2ImagesLeftRightOrient({imgs{1,1}; imgs{2,1}});

figure(3);
implottiling({img1; img2});
per

%%
[img1, img2, per] = align2ImagesLeftRightOrient({imgs{1,1}, imgs{1,2}});

figure(3);
implottiling({img1; img2});
per



%% Align 2 Images On Grid


%% RMS
aimgs ={imgs{1,1}; imgs{2,1}};
sh = align2ImagesOnGridByRMS(aimgs, 'overlap.max', 160, 'shift.max', 70)
sh = {[0,0]; sh};

figure(3); clf
plotAlignedImages(aimgs, sh);


%%

aimgs ={imgs{1,1}, imgs{1,2}};
sh = align2ImagesOnGridByRMS(aimgs, 'overlap.max', 160, 'shift.max', 70);
sh = {[0,0], sh};

figure(3); clf
plotAlignedImages(aimgs, sh)


%% Optimization


%% 

aimgs ={imgs{1,1}; imgs{2,1}};
sh = align2ImagesOnGridByOptimization(aimgs, 'overlap.max', 40, 'shift.max', 20)
sh = {[0,0]; sh};

figure(3); clf
plotAlignedImages(aimgs, sh)



%%

aimgs ={imgs{1,1}, imgs{1,2}};
sh = align2ImagesOnGridByOptimization(aimgs, 'overlap.max', 40, 'shift.max', 20);
sh = {[0,0], sh};

figure(3); clf
plotAlignedImages(aimgs, sh)



%% Hugin 

aimgs ={imgs{1,1}; imgs{2,1}};
sh = align2ImagesOnGridByHugin(aimgs, 'overlap.max', 50, 'shift.max', 30, 'project.filename', 'test.pto', 'project.cleanup', false, 'image.cleanup', true, 'image.filename', 'test')
sh = {[0,0]; sh};

figure(3); clf
plotAlignedImages(aimgs, sh)


%% Correlation


aimgs ={imgs{1,1}; imgs{2,1}};
sh = align2ImagesOnGridByCorrelation(aimgs, 'overlap.max', 50, 'shift.max', 30, 'project.filename', 'test.pto', 'project.cleanup', false, 'image.cleanup', false, 'image.filename', 'test')
sh = {[0,0]; sh};

figure(3); clf
plotAlignedImages(aimgs, sh)






%% AlignOnGrid - 2x2 grid

%% Generate 2x2 Grid

img = loadTestImage();
img = mat2gray(img);

size(img)

figure(1);
implot(img);

clear imgs
imgs{1,1} = img(1:230,1:250);
imgs{2,1} = img(200:end,20:260);
imgs{1,2} = img(10:240,230:end);
imgs{2,2} = img(230:end,240:end);

figure(2);
implottiling(imgs, 'link', false);



%% Global

sh = alignImagesOnGrid(imgs, 'overlap.max', 100, 'overlap.min', 1, 'shift.max', 50);
var2char(sh)

%figure(2);
%implottiling(imgs, 'link', false);

figure(3)
plotAlignedImages(imgs, sh)


%% Sequential 

sh = alignImagesOnGrid(imgs, 'overlap.max', 100, 'shift.max', 100, 'align', 'RMS', 'method', 'sequential');
var2char(sh)

%figure(2);
%implottiling(imgs, 'link', false);

figure(3)
plotAlignedImages(imgs, sh)


%% Primary

clc
sh = alignImagesOnGrid(imgs, 'overlap.max', 100, 'shift.max', 20, 'align', 'RMS', 'method', 'primary');
var2char(sh)

%figure(2);
%implottiling(imgs, 'link', false);

figure(3)
plotAlignedImages(imgs, sh)



%% AlignOnGird

%% 3x3 Tiling

img = loadTestImage();
img = mat2gray(img);

size(img)

figure(1);
implot(img);

imgs{1,1} = img(1:180,1:190);
imgs{2,1} = img(160:320,5:170);
imgs{3,1} = img(310:end,3:180);
imgs{1,2} = img(3:184,166:330);
imgs{2,2} = img(163:330,150:325);
imgs{3,2} = img(305:end,154:335);
imgs{1,3} = img(2:190, 320:end);
imgs{2,3} = img(175:335, 315:end);
imgs{3,3} = img(312:end, 312:end);


figure(2);
implottiling(imgs, 'link', false);


%% Global 

tic
sh = alignImagesOnGrid(imgs, 'overlap.max', 80, 'shift.max', 30, 'align', 'RMS');
toc
var2char(sh)

%figure(2);
%implottiling(imgs, 'link', false);

figure(3)
plotAlignedImages(imgs, sh)


%% Sequential 

tic
sh = alignImagesOnGrid(imgs, 'overlap.max', 80, 'shift.max', 30, 'align', 'RMS', 'method', 'sequential', 'anchor', [2,2]);
toc
var2char(sh)

%figure(2);
%implottiling(imgs, 'link', false);

figure(3)
plotAlignedImages(imgs, sh)


%% Primary

clc
tic
sh = alignImagesOnGrid(imgs, 'overlap.max', 100, 'shift.max', 50, 'align', 'RMS', 'method', 'primary');
toc
var2char(sh)

%figure(2);
%implottiling(imgs, 'link', false);

figure(3)
plotAlignedImages(imgs, sh)






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Images with arbitrary Overlapps


img = loadTestImage();
img = mat2gray(img);

size(img)

figure(1);
implot(img);

clear imgs
imgs{1} = img(1:280,1:290);
imgs{2} = img(160:320,5:170);
imgs{3} = img(200:400,3:380);
%imgs{4} = img(250:end,306:end);
%imgs{5} = img(163:330,150:325);

figure(2);
implottiling(imgs', 'link', false)

%% Full automatic alignment using Hugin

sh = alignImages(imgs, 'method', 'full')

figure(3);
plotAlignedImages(imgs, sh)


%% Overlapping Pairs 

clc
clear pairs
pairs(1).from = 1; pairs(1).to = 3;
pairs(2).from = 2; pairs(2).to = 3;

tic
sh = alignImages(imgs, 'pairs', pairs, 'overlap.max', 500, 'shift.max', 200);
toc
var2char(sh)

figure(3); clf
plotAlignedImages(imgs, sh)








