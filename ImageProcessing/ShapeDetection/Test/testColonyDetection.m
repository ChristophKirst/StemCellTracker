%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test Colony Detection Methods %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clear classes
close all
clc

initialize
bfinitialize



%% get data
texpr = tagexpr('./Test/Images/hESCells_Tiling/*.tif', 'tagnames', {'tile'})
is = ImageSourceTagged(texpr);

ist = ImageSourceTiled(is, 'tileshape', [4,4], 'tileformat', 'uy');

tl = ist.tiles;
size(tl)

figure(1); clf;
implottiling(tl)


% align

tic
ist.align('overlap.max', 120,  'overlap.min', 80, 'shift.max', 20)
toc

figure(1); clf
ist.plotAlignedImages


% stitch

img = ist.stitch('method', 'Mean');

figure(2); clf
implot(img);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Colony detection 


imgr = imresize(img, 0.5);
size(imgr)

imgrf = gaussianFilter(imgr, 20);
imgro = imclose(imgrf, strel('disk', 20));

figure(4); clf;
implottiling(cellfunc(@mat2gray, {imgr; imgrf; imgro}));


%% find threshold via histogram

figure(1); clf
hist(mat2gray(imgro(:)), 40)

th = thresholdFirstMin(mat2gray(imgro))


%%
imgm = mat2gray(imgro) > th;
implottiling(cellfunc(@mat2gray, {imgr; imgrf; imgro; imgm}));



%% identify regions and fit circle to them


imglab = bwlabeln(imgm);
imgsurf = impixelsurface(imglab);

%%
implottiling(cellfunc(@mat2gray, {imgr, imgrf, imgro; imgm, imglab, imgsurf}'));


%%
clc
stats= imstatistics(imglab, {'Volume', 'SurfacePixelIdxList'})

[stats.Volume]

%% select colonies with volume above certain rad
ids = find([stats.Volume] > pi*200^2)
ids = ids(1)

st = stats(ids)

st.BoundingBox
st.SurfacePixelIdxList;

is = zeros(size(imglab));
is(st.SurfacePixelIdxList) = 1;

figure(6); clf; 
implot(is)

%%

xy = imind2sub(size(imglab), st.SurfacePixelIdxList);

[c,r] = fitCircle(xy')

% add some safty margin:
r = r + 0.05 * r;

viscircles(c', r)



%%



%% find colonies

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


