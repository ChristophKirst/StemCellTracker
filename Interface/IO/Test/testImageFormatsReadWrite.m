%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test Image Formats IO changes %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clear classes
clc
close all

% note: we work from file: 
%    for a image either use imread, imshow, imwrite combination (yXZCT format)
%    or imreadBF, implot, imwriteBF qith XYCTZ format
%    -> images identical in orientation in matlab or exteranl fileviewer !
%    tested for tif !


%%
img = rand(5,7);
imgu = imrescale(img, 'class', 'uint8');

figure(1)
implottiling({img; imgu})


%% Export Import

% note: open external image viewer to check consistency woth orientation

%figure(1); clf
%imshow(img)

imwrite(img, 'test.tif')
imwrite(imgu, 'testu.tif');


imgir = imread('test.tif')
imguir = imread('testu.tif')

size(img)
size(imgir)
size(imguir)


figure(1); clf
imsubplot(1,3,1)
imshow(imgu)
imsubplot(1,3,2)
imshow(imgir)
imsubplot(1,3,3)
imshow(imguir)


%% info

iir = imfinfo('test.tif')
iuir = imfinfo('testu.tif')

%% 
clc
size(img)

size(impqlpermute(img, 'pqlct' , 'yplct'))


figure(1); clf
implot(img);

imwriteBF(img, 'test.tif')
imwriteBF(imgu, 'testu.tif');

imgbf = imreadBF('test.tif')
imgubf = imreadBF('testu.tif')

figure(1); clf
implottiling({imgu; imgbf; imgubf}, 'link', false)

%%
size(img)

{class(img), class(imgu)}
{class(imgir), class(imguir)}
{class(imgbf), class(imgubf)}


%% info
iibf = imwriteBFInfo('test.tif')
iiubf = imwriteBFInfo('testu.tif')


%%

delete('test.tif', 'testu.tif')



