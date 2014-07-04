%%%%%%%%%%%%%%%%%%%%%
%%% Test Stiching %%%
%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
initialize

%%

%%%%%%%%%%%%%%%%%%%%%
%%% 2D - 2 Images %%%
%%%%%%%%%%%%%%%%%%%%%

%% Load and Split a Test Image

img = imread('./Test/Images/hESCells_DAPI.tif');
img = mat2gray(img);
size(img)


%horizontal
img1 = img(1:303, 1:end-30);
img2 = img(260:end, 18:end);

%horizontal
img1 = img(1:303, 1:end-2);
img2 = img(260:end, 3:end);

%vertical
%img1 = img(10:end, 1:303);
%img2 = img(1:end-10, 260:end);


figure(1); clf
implot(img)

figure(2); clf
imsubplot(1,2,1)
implot(img1);
imsubplot(1,2,2);
implot(img2);

alignfigures


%% Load two images

img1 = imread('./Test/Images/hESCells_Tiling/W1F034T0001Z05C1.tif');
img2 = imread('./Test/Images/hESCells_Tiling/W1F035T0001Z05C1.tif');

img1 = impqlpermute(img1, 'yx', 'pq');
img2 = impqlpermute(img2, 'yx', 'pq');

size(img1)

figure(1); clf
imsubplot(1,2,1);
implot(double(img1));
imsubplot(1,2,2);
implot(img2);


%% align2ImagesByOptimization

% this routine is more robust to perturbations in seconadary directions
% choose overlap close to and above the expected average

tic 
sh = align2ImagesByOptimization(mat2gray(img1), mat2gray(img2), setParameter('overlap.max', 50, 'direction', 'lr'))
%sh = align2ImagesByOptimization(mat2gray(img1), mat2gray(img2), setParameter('overlap.max', 130, 'direction', 'lr'))
toc

figure(3); clf;
plot2AlignedImages(mat2gray(img1), mat2gray(img2), sh)


%% align2ImagesBySequentialShift 

% this routine is faster for good aligned images

tic
%sh = align2ImagesBySequentialShifts(mat2gray(img1), mat2gray(img2), setParameter('overlap',[50, 10, 10], 'minoverlap', [30, -10, -10], 'direction', 'lr'))
sh = align2ImagesBySequentialShifts(mat2gray(img1), mat2gray(img2), setParameter('overlap.max', 100, 'overlap.min', 30, 'shift.max', 10, 'shift.min',-10, 'direction', 'lr'))
toc

figure(3); clf;
plot2AlignedImages(mat2gray(img1), mat2gray(img2), sh)


%% stitch2ImagesByOverwrite


imgst = stitch2ImagesByOverwrite(img1, img2, sh);

figure(4); clf;
implot(imgst)


%% compare stitch2ImagesByXXX

imgst = stitch2ImagesByOverwrite(img1, img2, sh);
imgst2 = stitch2ImagesByMean(img1, img2, sh);
imgst3 = stitch2ImagesByMin(img1, img2, sh);
imgst4 = stitch2ImagesByMax(img1, img2, sh);

figure(4); clf;
implottiling({mat2gray(imgst), mat2gray(imgst2); mat2gray(imgst3),  mat2gray(imgst4)})



%% test alignment in z direction:

img1 = imread('./Test/Images/hESCells_Stack/W1F127T0001Z04C1.tif');
img2 = imread('./Test/Images/hESCells_Stack/W1F127T0001Z05C1.tif');

figure(1); clf;
implottiling({img1, img2});


sh = align2ImagesByOptimization(img1, img2, setParameter('overlap.max', size(img1, 1)))

figure(2); clf
plot2AlignedImages(mat2gray(img1), mat2gray(img2), sh)


%%
sh = align2ImagesBySequentialShifts(img1, img2, setParameter('overlap.max', size(img1, 1), 'overlap.min', size(img1,1)-5, 'shift.max', 5, 'shift.min', -5))

figure(3); clf
plot2AlignedImages(mat2gray(img1), mat2gray(img2), sh)



%% Load two images in y direction

img1 = imread('./Test/Images/hESCells_Tiling/W1F034T0001Z05C1.tif');
img2 = imread('./Test/Images/hESCells_Tiling/W1F038T0001Z05C1.tif');

img1 = impqlpermute(img1, 'yx', 'pq');
img2 = impqlpermute(img2, 'yx', 'pq');

size(img1)

figure(1); clf
imsubplot(2,1,1);
implot(img1);
imsubplot(2,1,2);
implot(img2);

sh = align2ImagesBySequentialShifts(img1, img2, setParameter('overlap.max', 140, 'overlap.min', 90, 'direction', 'tb'))

figure(3); clf
plot2AlignedImages(mat2gray(img1), mat2gray(img2), sh)



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2D - Grid of Images %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load linear grid of images

n1 = 4;
n2 = 4;
clear imgs

for i = 1:n1
   for j = 1:n2
      pos = 33 + (i-1) + 4 * (j-1);
      fn = ['./Test/Images/hESCells_Tiling/W1F' num2str0(pos,3), 'T0001Z05C1.tif'];
      imgs{i,j} = imread(fn);
      imgs{i,j} = impqlpermute(imgs{i,j}, 'yx', 'pq');
   end
end
imgs = imgs';
   
figure(1); 
clf
implottiling(imgs)

%%
imgsrow = imgs(2, :)

figure(2); clf
implottiling(imgsrow)


%% align image gird

shifts = alignImages(imgs, setParameter('method' , 'Optimization', 'overlap.max', 120))


%% align image gird

shifts = alignImages(imgsrow, setParameter('method' , 'SequentialShifts', 'overlap.max', 120, 'overlap.min', 90))

%% plot

clf
plotAlignedImages(imgsrow, shifts)


%% stitch 

img = stitchImagesByMin(imgsrow, shifts);

figure(1); clf; imcolormap('gray')
implot(img)


%%

shifts = alignImages(imgs, setParameter('method' , 'SequentialShifts', 'overlap.max', 120, 'overlap.min', 90))

figure(7)
clf
plotAlignedImages(imgs, shifts)


%% stitchImages

img = stitchImagesByOverwrite(imgs, shifts);

figure(8)
implot(img/2^16)


%% compare stitching methods

img1 = stitchImages(imgs, shifts, struct('method', 'Min'));
img2 = stitchImages(imgs, shifts, struct('method', 'Max'));
img3 = stitchImages(imgs, shifts, struct('method', 'Mean'));
img4 = stitchImages(imgs, shifts, struct('method', 'Overwrite'));


figure(10)
implottiling({img1/2^16, img2/2^16; img3/2^16, img4/2^16})





