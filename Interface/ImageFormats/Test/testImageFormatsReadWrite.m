%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test Image Formats IO changes %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clear classes
clc
close all

%%
img = rand(5,7);
imgu = imrescale(img, 'class', 'uint16');

figure(1)
implottiling({img; imgu})


%% Export Import


imwrite(img, 'test.tif')
imwrite(imgu, 'testu.tif');


%%

imgir = imread('test.tif')
imguir = imread('testu.tif')

imgbf = imread_bf('test.tif')
imgubf = imread_bf('testu.tif')



%% info

iir = imfinfo('test.tif')
iuir = imfinfo('testu.tif')

iibf = imread_bf_info('test.tif')
iiubf = imread_bf_info('testu.tif')




