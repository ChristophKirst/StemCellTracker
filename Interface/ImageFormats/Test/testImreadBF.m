%% example


clc
clear all
close all
initialize
bfinitialize

%%

info = imread_bf_info('./Test/Images/hESCells_tripple_reporter.TIF')


%%

dat = imread_bf('./Test/Images/hESCells_DAPI.tif');

figure(1)
implot(mat2gray(dat));

min(dat(:))
max(dat(:))

unique(dat(:));


%%
clc
info = imread_bf_info('./Test/Images/color.jpg')

%%

dat = mat2gray(imread_bf('./Test/Images/color.jpg', 'squeeze', false));
size(dat)

%%

imformat(dat)


%%
clc
datp = impqlpermute(dat, 'yxc', 'pqc');

figure(1)
implot(datp);




%%
clc
info = imread_bf_info('./Test/Images/hESCells_100x.tif')

%%

dat = mat2gray(imread_bf('./Test/Images/hESCells_100x.tif', 'squeeze', false));
size(dat)

%%
clc
datp = impqlpermute(dat, 'yxc', 'pqc');

figure(1); clf
implot(datp);








%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reading of various Formats

fname = './Test/Images/hESCells_Colony.zvi';

info = imread_bf_info(fname)

info(2)

%%

info = imread_bf_info(fname)
info(4)

%%
data1a = imread_bf(fname, 'series', 1, 'channel', 1);
size(data1a)

data1b = imread_bf(fname, 'series', 2, 'channel', 1);
data1c = imread_bf(fname, 'series', 3, 'channel', 1);
data1d = imread_bf(fname, 'series', 4, 'channel', 1);

%%
figure(5); clf;
implottiling({data1a, data1b; data1c, data1d})


%% 
data2 = imread_bf(fname);

size(data2)

%%
figure(5); clf;
implottiling({data1a; data2(:,:,1)})





