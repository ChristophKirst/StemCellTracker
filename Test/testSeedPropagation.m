%% Test Seed Propagation

clc
close all
clear all


%% test 2D 

img = zeros(100,100);
disk = mat2gray(fspecial2('disk', [20, 50], 0, 1));

seed = mat2gray(fspecial2('disk', 4));
imgmax =imreplace( imreplace(img, seed, [35, 45]), seed, [40, 80]);
img =imreplace(img, disk, [30, 40]);
img = imreplace(img, disk, [50,10]);


figure(1)
clf
imshow(mat2gray(img + imgmax))

figure(2)
imshow(img);


%% test propagation 
param.lambda = 0.03;
param.cutoff.difference = 0.5;
param.avaraging.ksize = 1;

[prop, dist] = seedPropagation(img, imgmax, ones(size(img)) > 0, param);

figure(3)
imsubplot(1,2,1);
imshow(prop)
imsubplot(1,2,2);
imshow(dist)


%%

clc
clear all
close all


%% test 3D

img = zeros(100,110,30);
disk = mat2gray(fspecial3('disk', [20, 50, 10]));


seed = mat2gray(fspecial3('disk', 4));
imgmax =imreplace( imreplace(img, seed, [35, 45, 11]), seed, [40, 75, 12]);

img =imreplace(img, disk, [30, 40, 7]);


% img = zeros(10,8,9);
% disk = fspecial3('disk', [5, 5, 5]) > 0;
% 
% sd = fspecial3('disk', 2) >0;
% imgmax = imreplace(img, sd, [3, 4, 5]);
% img =imreplace(img, disk, [2, 2, 3]);


figure(1)
clf
imsubplot(1,2,1);
immontage(img + 2* imgmax);
%implottiling(mat2gray(img + imgmax))

imsubplot(1,2,2);
immontage(img);



clc
unique(img + 2* imgmax)'
max(img(:))'
min(img(:))'
unique(img(:))'


%% test propagation

param.lambda = 0.05;
param.cutoff.difference = 0.5;
param.averaging.ksize = 3;

[prop, dist] = seedPropagation(img, imgmax, param);

figure(2)
clf
imsubplot(1,2,1);
immontage(prop);
imsubplot(1,2,2);
immontage(mat2gray(dist));

