%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test Hugin Interface %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
close all
clear all

initialize
hiinitialize
bfinitialize

%%
img = loadTestImage();

img1 = img(1:300, 1:300);
img2 = img(200:end, 1:300);
img3 = img(1:300, 200:end);
img4 = img(200:end, 200:end);

imgs= {img3, img4; img1, img2};
imgsg = cellfun(@mat2gray, imgs, 'UniformOutput', false);

figure(1); clf
implottiling(imgs)



%% Align Images using RMS method as reference

sh = align2ImagesByRMS(imgs{1}, imgs{2}, 'overlap.max', 500);

figure(2); clf;
plot2AlignedImages(imgs{1}, imgs{2}, sh)


%% Align via Hugin

sh = hialign(imgs, 'project.filename', 'test.pto', 'project.read', true, 'project.cleanup', false)
%sh = hipto2shifts(sh);
sh{2}

figure(3); clf;
plotAlignedImages(imgs, sh)


%% Stitch Via Hugin

imgst = histitch(imgs, sh);

figure(4); clf;
implot(imgst);


%% Align Real Data via Hugin

for t = 1:4
   imgs{t} = imread_bf('./Test/Images/hESCells_Colony.zvi', 'series', t, 'channel', 1);
   imgs{t} = imgs{t} - imopen(imgs{t}, strel('disk', 75));
end
imgs = {imgs{3}, imgs{4}; imgs{1}, imgs{2}};

figure(1); clf;
implottiling(imgs)

%%
sh = hialign(imgs, 'project.filename', 'test', 'project.cleanup', true, 'image.cleanup', true, 'image.filename', 'test')

figure(2); clf;
plotAlignedImages(cellfunc(@mat2gray, imgs), sh)


%% Generation of pto files

ptofile = hiptogenerate(imgs, 'project.filename', 'test.pto', 'project.read', false)
pto = hiparsepto(ptofile)

%%
clc 

sh = hipto2shifts(ptofile);
var2char(sh)

%%
isizes = cellfunc(@size, imgs);
pto2 = hishifts2pto(sh, isizes);

var2char({[pto.v], [pto.TrX], [pto.TrY]})
var2char({[pto.v], [pto2.TrX], [pto2.TrY]})



%% Align Images and only Stitch with Hugin

img = loadTestImage();

img1 = img(1:300, 1:300);
img2 = img(200:end, 1:300);
img3 = img(1:300, 200:end);
img4 = img(200:end, 200:end);

imgs= {img1, img3; img2, img4};
imgsg = cellfunc(@mat2gray, imgs);

figure(1); clf
implottiling(imgs, 'link', false)


sh = alignImagesOnGrid(imgs, 'overlap.max', 120, 'overlap.min', 80)


figure(1); clf
plotAlignedImages(imgs, sh)


%%

st = histitch(imgs, sh, 'temporary.cleanup', false);

figure(2); clf
implot(st)

%%

i1 = imread_bf('/tmp/tp06f22ffe_1917_4c3f_a26d_01a8bc2bd5fa0001.tif');
i2 = imread_bf('/tmp/tp06f22ffe_1917_4c3f_a26d_01a8bc2bd5fa0002.tif');

%st = stitchImages({i1, i2}, sh(1:2))

figure(1)
implottiling(cellfunc(@mat2gray,{i1(:,:,1); i2(:,:,1)}), 'link', false)


figure(2)
implottiling(cellfunc(@mat2gray,{i1(:,:,2); i2(:,:,2)}), 'link', false)



%%

i1 =imread_bf_info('/tmp/tp318a4c12_f82b_429c_9195_5ebcd5a8c8e60001.tif');
i2 =imread_bf_info('/tmp/tp318a4c12_f82b_429c_9195_5ebcd5a8c8e60002.tif');


i1.imetadata

i2.imetadata

%%

st2 = stitchImages(imgs, sh, 'method', 'Mean');

figure(3); clf
implottiling({st; st2}, 'link', false)


%%








%% Photometric optimization

it = imread('test0001.tif');
size(it)



%% generate project


[pto, ifiles] = hiptogenerate(imgs, sh, 'project.filename', 'test.pto', 'image.filename', 'test')


%% optimize 


system('autooptimiser -m -o popt.pto test.pto')




%% Photometric Correction of aligned Images


%% Generate 2x2 grid with intensity fluctuations 

img = loadTestImage();
img = mat2gray(img);

size(img)

figure(1);
implot(img)

clear imgs
imgs{1,1} = img(1:230,1:250);
imgs{2,1} = 0.5 * img(200:end,20:260);
imgs{1,2} = img(10:240,230:end);
imgs{2,2} = 1.5 * img(230:end,240:end);

figure(2); imcolormap('gray')
implottiling(imgs, 'link', false);






