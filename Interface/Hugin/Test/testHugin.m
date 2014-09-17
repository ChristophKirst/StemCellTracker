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

%% check

% i1 = imread_bf('/tmp/tp06f22ffe_1917_4c3f_a26d_01a8bc2bd5fa0001.tif');
% i2 = imread_bf('/tmp/tp06f22ffe_1917_4c3f_a26d_01a8bc2bd5fa0002.tif');
% 
% i1.imetadata
% i2.imetadata
% 
% 
% st = stitchImages({i1, i2}, sh(1:2))
% 
% figure(1)
% implottiling(cellfunc(@mat2gray,{i1(:,:,1); i2(:,:,1)}), 'link', false)





%% Photometric optimization

it = imread('test0001.tif');
size(it)



%% generate project


[pto, ifiles] = hiptogenerate(imgs, sh, 'project.filename', 'test.pto', 'image.filename', 'test')


%% optimize 


system('autooptimiser -m -o popt.pto test.pto')




%% Photometric Correction of aligned Images

clc
close all
clear all

initialize
hiinitialize
bfinitialize


%% Generate 2x2 grid with intensity fluctuations 

img = loadTestImage();
img = mat2gray(img);

size(img)

figure(1);
implot(img)

clear imgs
imgs{1,1} = img(1:230,1:250);
imgs{2,1} = 0.95 * img(200:end,20:260);
imgs{1,2} = img(10:240,230:end);
imgs{2,2} = 1.05 * img(230:end,240:end);

figure(2); imcolormap('gray')
implottiling(imgs, 'link', false);


%% align

sh = alignImagesOnGrid(imgs, 'overlap.max', 50, 'overlap.min', 1, 'shift.max', 50, 'alignment', 'RMS')

figure(1); clf
plotAlignedImages(imgs, sh)


%% generate pto project 


[pto, ifiles] = hiptogenerate(imgs, sh, 'project.filename', 'test.pto', 'image.filename', 'test')

%%





%% 

v Ra0 Rb0 Rc0 Rd0 Re0 Vb0 Vc0 Vd0
v Eev1 Eev2 Eev3 Eev4


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - real data


te= tagexpr('./Test/Images/hESCells_GCamp_Vignetting/*.tif', 'tagnames', {'field'})

is = ImageSourceTagged(te);
is.setTagRange('field', {1,2,3,4,5,6});

ist = ImageSourceTiled(is, 'tileshape', [3,2], 'tileformat', 'uy');
ist.print

imgs = ist.tiles;

figure(1); clf
implottiling(imgs)


%%

ist.align

figure(2);
ist.plotAlignedImages


%%
%sh = hialign(imgs,'project.filename', 'test.pto', 'image.filename', 'test',  'project.cleanup', false, 'image.cleanup', false)

var2char(sh)
var2char(ist.imageShifts)

sh2 = ist.imageShifts;


%%

st = histitch(imgs, sh2);


%%
figure(6); clf
implot(imrescale(st, 'data.max', 20000))


%%

%sh = ist.imageShifts;
[pto, ifiles] = hiptogenerate(imgs, sh, 'project.filename', 'test.pto', 'image.filename', 'test')


%%

%system('pto_var --opt Ra0,Rb0,Rc0,Rd0,Re0,Va0,Vb0,Vc0,Vd0,Vx0,Vy0,Eev -o photo.pto test.pto')
system('pto_var --opt Va,Vb,Vc,Vd,Vx,Vy,Eev -o photo.pto test.pto')


%%

%imgsr = imrescaleall(imgs, 'class', 'uint16');
imgsr = imgs;

for i = 1:6
   img = imgsr{i};
   %imgc = imgray2color(img);
   
   
   mxval = immaxvalue(class(img));

   %imgo = cat(3, img, mxval * ones(size(img), 'like', img));

   imwrite_tiff(mat2gray(cat(3,img,img, img)), ['test', num2str0(i,4), '.tif'])
   %imwrite_tiff(cat(3,mat2gray(img),ones(size(img))), ['test', num2str0(i,4), '.tif'])
   
   %imwrite(imgs{i}, ['test', num2str0(i,4), '.tif'])
   %imwrite(imgs{i}, ['test', num2str0(i,4), '.tif'], 'WriteMode', 'append')
end

%% 

system('vig_optimize -o opt.pto photo.pto')


%%

pto = hiparsepto('opt.pto')
pto(1)

%%
ii = 2;

vi = hivignetting(size(imgs{ii}), pto(ii));

vi = vi + 0.00;


figure(7); clf
implot(vi);

max(vi(:))
min(vi(:))

figure(8); clf; 
implottiling(cellfunc(@mat2gray, {imgs{ii}; mat2gray(vi); imgs{ii}./vi}))




%%

system('autooptimiser -m -o auto.pto photo.pto')

%%

system('nona -z LZW -r ldr -m TIFF_m -o nona auto.pto')

%%
clc

for i = 1:4
   %imgg{i} = imread(['test0001 - test0004000', num2str(i-1), '.tif']);
   imgg{i} = imread(['nona000', num2str(i-1), '.tif']);
   size(imgg{i})
end
   
figure(7); clf;
implottiling(cellfunc(@(x) x(:,:,1), imgg)', 'link', false)

%%

imfinfo('test0001.tif')

%%


clc
for i = 1:4

   img = mat2gray(imgs{i});

t = Tiff(['test000', num2str(i), '.tif'], 'w'); 
tagstruct.ImageLength = size(img, 1); 
tagstruct.ImageWidth = size(img, 2); 
tagstruct.Compression = Tiff.Compression.None; 
%tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP; 
tagstruct.Photometric = Tiff.Photometric.MinIsBlack; 
tagstruct.BitsPerSample = 16; 
tagstruct.SamplesPerPixel = 1; 
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky; 
t.setTag(tagstruct); 
t.write(uint16(img));
%t.write(single(cat(3, img, img, img))); 
t.close();

end


%%

imfinfo('test0001.tif')

%%
in = imread_bf_info('test0001.tif')
in.imetadata

%%


img = imread('nona.tif');
size(img)


%%

for i = 1:4
   %imgg{i} = imread(['test0001 - test0004000', num2str(i-1), '.tif']);
   imgg{i} = imread(['nona000', num2str(i-1), '.tif']);
   size(imgg{i})
end
   

%%
figure(7); clf;
implottiling(cellfunc(@(x) x(:,:,1), imgg))

   
   



