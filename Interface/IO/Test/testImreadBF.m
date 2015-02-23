%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test imreadBF Bioformats %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all
initialize
bfinitialize

%%

[info, ir] = imreadBFInfo('./Test/Images/hESCells_tripple_reporter.TIF')

%%
clc

dat = imreadBF('./Test/Images/hESCells_DAPI.tif');

figure(1)
implot(mat2gray(dat));

min(dat(:))
max(dat(:))

unique(dat(:));


%%
clc
info = imreadBFInfo('./Test/Images/color.jpg')

%%

dat = mat2gray(imreadBF('./Test/Images/color.jpg', 'squeeze', true));
size(dat)

%%

imfrmtFormat(dat)


%%
clc
datp = imfrmtReformat(dat, 'YXC', 'XYC');

figure(1)
implot(datp);




%%
clc
info = imreadBFInfo('./Test/Images/hESCells_100x.tif')

%%

dat = mat2gray(imreadBF('./Test/Images/hESCells_100x.tif', 'squeeze', false));
size(dat)

%%
clc
datp = imfrmtReformat(dat, 'yxc', 'XYC');

figure(1); clf
implot(datp);








%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reading of various Formats

fname = './Test/Images/hESCells_Colony.zvi';

info = imreadBFInfo(fname)


%%
data1a = imreadBF(fname, 'S', 1, 'C', 1);
data1b = imreadBF(fname, 'S', 2, 'C', 1);
data1c = imreadBF(fname, 'S', 3, 'C', 1);
data1d = imreadBF(fname, 'S', 4, 'C', 1);

imgs = {data1a, data1b; data1c, data1d};

%%
size(data1d)
class(data1d)

%%
figure(5); clf;
implottiling(imgs)

figure(6); clf;
imgsr = imfrmtReformat(imgs, 'UV', 'Vu');
implottiling(imgsr)

%% 
imgs = imreadBF(fname, 'C', 1);

size(imgs)
size(imgs{1})

%%
imgs = reshape(imgs, [2,2]);

figure(5); clf;
implottiling(imgs)

figure(6); clf;
imgsr = imfrmtReformat(imgs, 'UV', 'Uv');
implottiling(imgsr)



%% Using the Image Reader

[info, ir] = imreadBFInfo('./Test/Images/hESCells_Colony.zvi')

%%
img =imreadBF(ir, 'C', 1, 'S', 1);
size(img)


%% ranges and sizes


[si, frmt] = imreadBFSize('./Test/Images/hESCells_Colony.zvi')
%[si, frmt] = imreadBFSize('./Test/Images/hESCells_Colony.zvi', 'xy')


%%

r = imreadBFRange('./Test/Images/hESCells_Colony.zvi', 'format', 'XYZ', 'y', 1:3)


