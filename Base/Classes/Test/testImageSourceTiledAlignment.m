%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ImageSourceTiled

clear all
clear classes
close all
clc

initialize
bfinitialize

texpr = tagexpr('./Test/Images/hESCells_Tiling/*.tif', 'tagnames', {'tile'})
is = ImageSourceTagged(texpr);

ist = ImageSourceTiled(is, 'tileshape', [4,4], 'tileformat', 'uy');


%%

tl = ist.tiles;
size(tl)

figure(1); clf;
implottiling(tl)


%%

tic
ist.align('overlap.max', 120,  'overlap.min', 80, 'shift.max', 20)
toc

%%
figure(1); clf
ist.plotAlignedImages


%%

img = ist.stitch('method', 'Mean');

figure(2); clf
implot(img);



%%


imgr = imresize(img, 0.5);
size(imgr)

imgrf = gaussianFilter(imgr, 20);
imgro = imclose(imgrf, strel('disk', 20));



figure(4); clf;
implottiling(cellfunc(@mat2gray, {imgr; imgrf; imgro}));


%%

figure(1); clf
hist(mat2gray(imgro(:)), 256)

%%
imgm = mat2gray(imgro) > 0.1;
implottiling(cellfunc(@mat2gray, {imgr; imgrf; imgro; imgm}));


%%

[centers, radii, metric] = imfindcircles(imgro,  [200 400], 'Sensitivity', 0.97, 'Method', 'TwoStage')

figure(6); clf
%implot(imgro);

imshow(mat2gray(imgro))
viscircles(centers,radii);

[~, id] = max(metric);

id = metric > 0.15
viscircles(centers(id,:), radii(id), 'EdgeColor', 'b')


%% aryehs approach:
% find peaks method to detect possible nucelar positions  -> alpha vol to detec colony -> works well -> go for it now, need to do it anyway

% methods to find circles ???




imgd = mat2gray(imgr);

imgf = medianFilter(imgd, 5);
imgf = mat2gray(imgf);

imgp =imextendedmax(imgf, 0.01);


figure(11); clf; 
implottiling({255 * imgd; 255 * imgf; imoverlay(imgd, imgp)})


%%
[p,q] = find(imgp);

X= [p,q];

clf
 alphavol(X,20,1);


