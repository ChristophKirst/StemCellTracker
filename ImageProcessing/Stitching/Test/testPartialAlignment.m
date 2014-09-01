%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test alignment of Black Images %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clear classes
close all
clc
initialize

%% Overlap Statistics

img = loadTestImage();

img1 = img(1:230,:);
img2 = img(200:end, 20:end-10);

p = AlignmentPair(img1, img2);
p.align()

figure(1); clf
p.plot()

%%

[o1, o2] = p.overlap();

figure(2); clf
implottiling({o1; o2})

%%
clc
stats = overlapStatisticsImagePair(p, 'overlap.max', 100, 'aligned', true);

stats.from
stats.to

%% non -alinged statistics (default)
clc
stats = overlapStatisticsImagePair(p, 'overlap.max', 100);

stats.from
stats.to


%% one black image
clc
p = AlignmentPair(img1, 0.01 * rand(100,100));

stats = overlapStatisticsImagePair(p);
stats.from
stats.to


%% overlap quality

q = overlapQuality(stats, 'threshold.var', 10)



%% determine connected components 







