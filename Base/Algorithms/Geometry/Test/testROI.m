%%%%%%%%%%%%%%%%%%%%%%
%%% ROIs           %%%
%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
clear classes

initialize

%%

r = ROIRectangle([2, 5; 2, 3])

r.toArray

r.volume

r.npixel

m = r.mask([10,10]);
figure(1); clf
implot(m);

%%
r
