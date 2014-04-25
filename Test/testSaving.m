%%%%%%%%%%%%%%%
% Test Saving %
%%%%%%%%%%%%%%%


%% 
clear all
close all
clc

initialize


%% Generate some test Data

isize = [100, 100];

imglab = syntheticLabeledImage(isize, 3, 5);
img = rand(isize);
figure(1);
implottiling({imcolorize(imglab), imgray2color(img, 'r')})

objs = label2DataObjects(imglab, img)

imgl = objs.labeledImage();
any(any(imglab- imgl))

objs.dataFields
[objs.value('MedianIntensity')]

%% Saving

save('./Test/Data/testdata.mat', 'objs')


%%
clear all
load('./Test/Data/testdata.mat')


%%

% setup an experiment























