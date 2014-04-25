%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test ImageJ interface %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initialize

%%
clear all
close all
javarmdynamicclasspath

%%
ijcompile(1)


%%
ijpath

%%
ijinitialize


%%

img = syntheticLabeledImage([30, 50, 10], 20, 10);
figure(1); clf
colormap gray
implot3d(img)
freezecolormap

%%
ijplot3d(img)


%%

figure(2)
imgc = imcolorize(img);
implot3d(imgc)

%%
ijplot3d(imgc)
