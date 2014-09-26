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


%% ROIPolygon


p = 3 + [1,1; 1,6; 3,9; 8,12 ]'

roi = ROIPolygon(p)


bb= roi.boundingbox
bb.toPixelArray


r = roi.mask([20,20]);
figure(1); clf
implot(r)


%%
clc
d = rand(20,30);
de= roi.extractdata(d);

de2 = immask(d, roi.mask(size(d)));

figure(2); clf
implottiling({d; de; de2}, 'link', false)






