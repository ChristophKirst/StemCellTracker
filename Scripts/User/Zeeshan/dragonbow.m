%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Dragonbow Cluster Detection %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

initialize
bfinitialize


%% Image Data

ddir = '/home/ckirst/Science/Projects/StemCells/Experiment/Data/DragonBow';
fh  = FileHandler('basedirectory', ddir, 'datadirectory', '.', 'resultdirectory', '.');
fh.printInfo
% 
% {fh.dirDataDirectory.name}
% {fh.dirDataDirectory.bytes}

%%

is = ImageSourceBF(fh.dataFile('5.tif')) 

is.setCellDataFormat('XY', 'C')
is.printInfo

%%
clc
d = is.cell;

figure(1); clf
implottiling(d)

%% denoise

d = is.cell;

dbkg = cellfunc(@(x) imopen(x, strel('disk', 100)), d);
dc   = cellfunc(@(x,y) x-y, d, dbkg);
dc   = cellfunc(@(x) double(x) / max(x(:)), dc);

figure(1); clf
implottiling({dc{:}}')

%% color image

imgc = cat(3, dc{:});
imgc = imgc/max(imgc(:));
figure(2); clf
implot(imgc)

%figure(3); clf
%implot(sum(imgc,3))

%dcc = cellfunc(@(x) imopen(x, strel('disk', 5)), dc);
dcc = cellfunc(@(x) filterGaussian(x, 35), dc);
imgcc = cat(3, dcc{:});

figure(3); clf
implot(imgcc)



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
