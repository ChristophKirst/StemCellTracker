%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Dragonbow Cluster Detection %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

initialize
bfinitialize

verbose = true;

%% Image Data

ddir = '/home/ckirst/Science/Projects/StemCells/Experiment/DragonBow';
fh   = FileHandler('basedirectory', ddir, 'datadirectory', '.', 'resultdirectory', '.');
fh.printInfo
% 
% {fh.dirDataDirectory.name}
% {fh.dirDataDirectory.bytes}

%%

is = ImageSourceBF(fh.dataFile('DRAGONbow_tile_c1+2+3+4.tif'));
%is.setCellDataFormat('XY', 'C');
is.setRange('X', 1000:3000, 'Y', 1000:2000);
is.printInfo


%%
clc
img = is.data / 255;
class(img)
size(img)

figure(1); clf
implot(img)


%%
imgf = filterBM(img, 'sigma', 10, 'profile', 'lc');
figure(2); clf
implottiling({img, imgf})

imgfLarge = imgf;

%% subimage for testing

imgf =imgfLarge;
imgf = imgf(500:1100, 1:500, :);
imgr = img(500:1100, 1:500, :);
figure(3); clf
implot(imgf)

%% SLIC

imglab = segmentBySLIC(imgf, 'superpixel', 2500, 'sigma', 5);
imglabs = impixelsurface(imglab);

figure(2); clf 
%colormap parula
imgSLIC = imoverlaylabel(imgf, imglabs, false)
implot(imgSLIC)


%% Masking and Intensity

imgI = imgf(:,:,1) + imgf(:,:,2) + imgf(:,:,3);
imgI = mat2gray(imgI);

imgmask = imgI > 0.05;
%imgmask = imopen(imgmask, strel('disk', 3));
%imgmask = postProcessSegments(bwlabeln(imgmask), 'volume.min', 50) > 0;

if verbose
   %max(img(:))
   figure(21); clf;
   set(gcf, 'Name', ['Masking'])
   implottiling({imoverlaylabel(imgI, double(imgmask));  imgmask})
end



%% Postprocess
clc
max(imglab(:))

stats = imstatistics(imglab, {'MinIntensity'}, imgI);
figure(7); clf; 
hist([stats.MinIntensity], 256)
hist(imgI(:), 256)


%%

%imgpost = postProcessSegments(imglab, imgI, 'intensity.min', 0.05, 'volume.min', 15, 'fillholes', false);
imglabm = immask(imglab, imgmask);
%imglabm = imlabelapplybw(imgpost, @(x) imopen(x, strel('disk', 2)));
imglabm = imlabelseparate(imglabm);
imglabm = imrelabel(imglabm);

max(imglab(:))
max(imglabm(:))

if verbose
   imgSp = impixelsurface(imglabm);
   figure(5); clf;
   imgNeurons = imoverlaylabel(imgf, imgSp, false)
   implot(imgNeurons);
end

%%

stats = imstatistics(imglabm, {'PixelIdxList', 'MinIntensity', 'MedianIntensity', 'Volume'}, imgI);

figure(8); clf; 
hist([stats.MinIntensity], 256)
hist(imgI(:), 256)

%%
ids = zeros(1,length(stats));
%ids = [stats.Extent] < 0.5; %* pi /4;
%ids = or(ids, [stats.FilledArea] > ([stats.Volume] + 10));
%ids = or(ids, 0.65 * ([stats.Perimeter]) / (2 * pi) >= sqrt([stats.Volume] / pi));
ids = or(ids, [stats.Volume] <= 20);
%ids = or(ids, [stats.MedianIntensity] < 0.01);
%ids = or(ids, [stats.MinIntensity] < 0.01);

%ids = ~ids;

imgpost = imglabm;
for i = find(ids)
   imgpost(stats(i).PixelIdxList) = 0;
end
imgpost = imrelabel(imgpost);
max(imgpost(:))


if verbose
   imgSp = impixelsurface(imgpost);
   figure(5); clf;
   implot(imoverlaylabel(imgf, imgSp, false));
end


%% Colorize

rp = imstatistics(imgpost, {'Volume', 'PixelIdxList'});

imgR = imgf(:,:,1);
imgG = imgf(:,:,2);
imgB = imgf(:,:,3);

imgLR = zeros(size(imgR));
imgLG = zeros(size(imgG));
imgLB = zeros(size(imgB));

l = zeros(length(rp), 3);

for i = 1:length(rp)
   px = rp(i).PixelIdxList;
   r = mean(imgR(px));
   g = mean(imgG(px));
   b = mean(imgB(px));
   
   imgLR(px) = r;
   imgLG(px) = g;
   imgLB(px) = b;
   
   l(i,:) = [r,g,b]; 
end


imgL = cat(3, imgLR, imgLG, imgLB);

figure(10); clf
implottiling({imgf; imgL})
   

%% spatial distances

nl = max(imgpost(:));
distL = cell(1,nl);
st = imstatistics(imgpost, {'PixelIdxList'});

for l = 1:nl
   distL{l} = bwdist(imgpost == l);
end

distS = zeros(nl);
for i = 1:nl
   for j = i+1:nl
      pxi = st(i).PixelIdxList;
      dL = distL{j};
      distS(i,j) = min(dL(pxi));
      distS(j,i) = distS(i,j);
   end
end


figure(15);
implot(distS)
   

%% Cell line detection via clustering

% Lab color space
lb = rgb2lab(l);


figure(7); clf
subplot(1,2,1)
plot3(l(:,1), l(:,2), l(:,3), 'k+')

subplot(1,2,2)
plot3(lb(:,1), lb(:,2), lb(:,3), 'k+')



%%
lpy = numpyFromMat(lb);

cl = py.scipy.cluster.vq.kmeans(lpy,10);
cl = numpyToMat(cl);

dli = distanceMatrix(lb', cl{1}');
[ml, pl] = min(dli, [], 2);

clc = lab2rgb(cl{1})

imgLR = zeros(size(imgR));
imgLG = zeros(size(imgG));
imgLB = zeros(size(imgB));

for i = 1:length(pl)
   px = rp(i).PixelIdxList;

   r = clc(pl(i), 1);
   g = clc(pl(i), 2);
   b = clc(pl(i), 3);
   
   imgLR(px) = r;
   imgLG(px) = g;
   imgLB(px) = b; 
end


imgL2 = cat(3, imgLR, imgLG, imgLB);
imgL2 = imgL2 / max(imgL2(:));

figure(15); clf
implottiling({imgf; imgL; imgL2})
 
 
 
 %% Save Images
 
imwriteBF(imgr, [ddir, '/DragonBow_Raw.jpg']);
imwriteBF(imgf, [ddir, '/DragonBow_Filtered.jpg']);
imwriteBF(imgSLIC, [ddir, '/DragonBow_SLIC.jpg']);
imwriteBF(imgNeurons, [ddir, '/DragonBow_Neurons.jpg']);
imwriteBF(imgL, [ddir, '/DragonBow_Segments.jpg']);
imwriteBF(imgL2, [ddir, '/DragonBow_ColorCluster.jpg']);


%% super pixel segmentation



 





%% color image

dt = cellfunc(@(x) imclose(x, strel('disk', 5)), dc);
imgt = cat(3, dt{3}, dt{1}, dt{2});
imgt = 1-imgt/max(imgt(:));
figure(3); clf
implot(imgt)



%%
%figure(3); clf
%implot(sum(imgc,3))

%dcc = cellfunc(@(x) imopen(x, strel('disk', 5)), dc);
dcc = cellfunc(@(x) filterDisk(x, 35), dc);
imgcc = cat(3, dcc{:});

figure(3); clf
implot(imgcc)


imgh = rgb2hsv(imgcc);
figure(6)
implottiling({imgh(:,:,1); imgh(:,:,2); imgh(:,:,3)})


figure(7)
h = imgh(:,:,1);
subplot(1,2,1)
hist(h(:), 200)

subplot(1,2,2)
h = imgh(:,:,3);
hist(h(:), 100)

%%
th = thresholdFirstMin(h(:), 'delta', 100)
msk = imgh(:,:,3) >  th;
msk  = imdilate(msk, strel('diamond', 1));
msk = msk > 0;

figure(9); clf
implottiling({imgc; msk})


figure(10); clf
hh = imgh(:,:,1);
hh = hh(msk);
hist(hh(:), 200)
[hd, hx] = hist(hh(:), 200);

%% 

tt = findTroughs(hd, 400);
tt = hx(tt)

%%
% define clusters

pos = [204, 402; 246, 890; 534, 500; 788, 754; 867, 338]';

for i = 1:size(pos,2);
   cols(i,:) = imgcc(pos(1,i), pos(2,i), :);
end

%%
% imgi = sum(imgcc,3);
% figure(8); clf;
% hist(imgi(:), 100)
% 
% th = thresholdFirstMin(imgi(:))
% 
% msk = imgi > th;
% figure(14); clf
% implottiling({imgc; msk})


%%
dccl = cellfunc(@(x) x(:), dcc);
iccl = cat(2, dccl{:});

di = distanceMatrix(iccl', cols');
[m, p] = min(di, [], 2);
%%

imglab = zeros(size(dcc{1}));
imglab(:) = p;
imglab = immask(imglab, msk);

figure(14); clf
colormap jet
implottiling({imgc; imglab})



ht = findTroughs(hh, 300)
tt = xx(ht)
tt = [tt, 0.865];
tt = [0, tt, max(dh(:))+eps];

length(tt)

%% labeled image

lab = zeros(size(dcc{1}));
dhh = imgh(:,:,1);

for i = 1:length(tt)-1
   lab(and(dhh >= tt(i), dhh < tt(i+1))) = i;
end
lab(~msk) = 0;

figure(15); clf
implottiling({imgc; lab})























%% rescales image

dr = cellfunc(@(x) imresize(x, 1), dc);
drl = cellfunc(@(x) x(:), dr);

imgr = cat(3, dr{:});

figure(4); clf
implot(imgr)


%% clutering in rgb space

dnl = cellfunc(@(x) x(:), dcc);
dnl = cat(2, dnl{:});
dpy = numpyFromMat(dnl);
%dpy.shape

cl = py.scipy.cluster.vq.kmeans(dpy,6);
cl = numpyToMat(cl);

%%
di = distanceMatrix(dnl', cl{1}');
[m, p] = min(di, [], 2);
%%

imglab = zeros(size(dr{1}));
imglab(:) = p;

figure(14); clf
colormap jet
implottiling({imgr; imglab})



%% plot data as intensities

figure(3)
plot3(d3{:}, '*')
xlabel('R'); ylabel('G'); zlabel('B')
figure(4); clf
dn = cellfunc(@(x) x/ max(x(:)), dr);
implot(cat(3, dn{:}))

figure(5); clf
implottiling(reshape([{cat(3, dn{:})}; dn], 2,2))



%% hsv

figure(12); clf
ndhist(dcc{1}(:), dcc{2}(:))



%% Lab color space

imgl = rgb2lab(imgr);

figure(8)
dl = cellfunc(@(x) imgcc(:,:,x), {1;2;3});
implottiling(dl)

%% background detection
figure(9); clf
hist(dl{1}(:), 50)

th = thresholdFirstMin(dl{1}, 'delta', 20)
%th = 24;

msk = dl{1} >  th;
msk = imclose(msk, strel('disk', 10));

figure(10); clf
implottiling({imgr;  dl{1}; msk * 255},'link', false)


%% clustering on Lab space

dll = cellfunc(@(x) x(:), dl(2:3));
dll = cat(2, dll{:});
dlpy = numpyFromMat(dll);
%dpy.shape

cll = py.scipy.cluster.vq.kmeans(dlpy,15);
cll = numpyToMat(cll);

dli = distanceMatrix(dll', cll{1}');
[ml, pl] = min(dli, [], 2);


imglab = zeros(size(dr{1}));
imglab(:) = pl;

figure(14); clf
colormap jet
implottiling({imgr; imglab})




%% Hue color space

imgl = rgb2hsv(imgc);

figure(8)
dl = cellfunc(@(x) imgcc(:,:,x), {1;2;3});
implottiling(dl)

%% background detection
figure(9); clf
hist(dl{1}(:), 50)

th = thresholdFirstMin(dl{1}, 'delta', 20)
%th = 24;

msk = dl{1} >  th;
msk = imclose(msk, strel('disk', 10));

figure(10); clf
implottiling({imgr;  dl{1}; msk * 255},'link', false)


%% clustering on Lab space

dll = cellfunc(@(x) x(:), dl(2:3));
dll = cat(2, dll{:});
dlpy = numpyFromMat(dll);
%dpy.shape

cll = py.scipy.cluster.vq.kmeans(dlpy,15);
cll = numpyToMat(cll);

dli = distanceMatrix(dll', cll{1}');
[ml, pl] = min(dli, [], 2);


imglab = zeros(size(dr{1}));
imglab(:) = pl;

figure(14); clf
colormap jet
implottiling({imgr; imglab})





%% Hue color space

dn = cat(3, dc{:});
dn = dn / max(dn(:));

dnc = rgb2hsv(dn);

figure(8)
dncc = cellfunc(@(x) dnc(:,:,x), {1;2;3});
implottiling(dncc)

%% background detection
figure(9); clf
hist(dncc{3}(:), 50)

th = thresholdFirstMin(dncc{3}, 'delta', 20)
%th = 24;

msk = dncc{3} >  th;
msk = imclose(msk, strel('disk', 10));

figure(10); clf
implottiling({dd/max(dd(:));  dncc{3}; msk * 255},'link', false)


dh = dncc{1};
dh = dh(dncc{3} > th);

figure(11); clf
hist(dh(:), 100)
[hh, xx] = hist(dh(:), 100);



%% 

% peaks antrhoughs on hue

ht = findTroughs(hh, 300)
tt = xx(ht)
tt = [tt, 0.865];
tt = [0, tt, max(dh(:))+eps];

length(tt)

%% labeled image

dhh = dncc{1};

lab = zeros(size(dncc{1}));
for i = 1:length(tt)-1
   lab(and(dhh >= tt(i), dhh < tt(i+1))) = i;
end
lab(dncc{3} < th) = 0;

figure(15); clf
implottiling({dn; lab})

%% k means clustering to detect patches
%(is ther a way to include spatial info ??)

del = dncc{1} <  th;
a = dncc{2};
b = dncc{3};

ab2 = a.*a + b.*b;

a(del) = [];
b(del) = [];

figure(12); clf
ndhist(a,b)




%%

dd1 = dn;
aa = dncc{2};
bb = dncc{3};
for i = 1:3
   q = dd1(:,:,i);
   q(~and(aa > 50, bb < -20)) = 0;
   dd1(:,:,i) = q;
end

figure(13); clf
implottiling({dd1; dd})
