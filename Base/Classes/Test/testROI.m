%%%%%%%%%%%%%%%%%%%%%%
%%% ROIs           %%%
%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
clear classes

initialize

%% ROIRectangle
clc
r = ROIRectangle([2, 3], [7, 5])

r.toArray
r.toPixelArray

%%
clc
r.volume
r.npixel

mm = r.mask([10,10]);
total(mm)

%%
m = r.mask([10,10]);
figure(1); clf
implot(m);

%% ROIDisk

clc
r = ROIDisk([15,20], 10)

m = r.mask([30,40]);
r.npixel
total(m)


figure(2); clf
implot(m)


%% ROIMask

clc
r = ROIMask(rand(20,30)>0.5, 10)

m = r.mask;
r.npixel
total(m)


figure(2); clf
implot(m)



