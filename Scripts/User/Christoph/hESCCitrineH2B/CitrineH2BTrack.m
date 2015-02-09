%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Segment and Track Citrine H2B Movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clear classes
close all
clc
initialize
bfinitialize
addpath('./Scripts/User/Christoph');
verbose = false;
initializeParallelProcessing(12) % number of processors
%% Initialize Files
fns = '/var/run/media/ckirst/ChristophsData/Eric/H2BCitrine/Analysis/ColonyData/ColonyW10T<T,4>Z<Z,2>C1.tif';
is = ImageSourceFiles(fns);
%is.setCellDataFormat('XYZ', 'T');
is.printInfo
%% 3D Segmentation
parfor time = 1:10
stackraw = is.cell('T', time);
stackraw = cat(3, stackraw{:});
%%
%
% if verbose
% figure(1)
% clf
% boxr = [1 1 1];
% set(gcf, 'Name', 'Raw Stack');
% implot3d(stackraw, 'BoxRatios', boxr, 'Texture', 'z');
% ijplot3d(stackraw, 'PixelDepth', 1)
% end
%figure(1); clf
%hist(stackraw(:), 2560)
th = 2200;
imgmask = stackraw > th;
imgmask = imclose(imgmask, strel('disk', 5));
if verbose
figure(2); clf
set(gcf, 'Name', ['Time: ', num2str(time), ' Mask']);
implottiling(imgmask, 'tiling', [4,3])
drawnow
end
imgseg = CitrineH2BSegment3D(stackraw, imgmask, verbose, time);
[objects, stats] = label2Objects(imgseg, stackraw, 'time', time);
frame(time) = Frame('objects', objects);
% save image data
imgs{time} = {stackraw, imgseg};
end
%% plot
maxtime = length(imgs);
movie = zeros([is.dataSize, is.cellSize('Z'), 3, maxtime]);
for t = 1:maxtime
imgovl = imoverlaylabel(imgs{t}{1}, impixelsurface(imgs{t}{2}), false);
%implottiling(imgovl, 'tiling', [4,3]);
movie(:,:,:,:,t) = imgovl;
end
size(movie)
%%
%ijplot5d(movie)
figure(12); clf
implottiling(movie(:,:,:,:,2), 'tiling', [4,3])
%% Tracking
param.optimize = true;
param.print.optimization = true;
[match, cost] = matchObjects(frame(1), frame(2), param);
figure(18)
clf
subplot(1,2,1)
plotMatchedObjects(match)
title('matches')
subplot(1,2,2)
plotMatchedCostMatrix(match, cost)
title('cost matrix')
%% match frames
%method
param.optimize = true; % optimize match using an additional optimal coordinate transformation and second matching
% cost matrix
param.cost.creation = []; % costs for creating new objects ([])
param.cost.deletion = []; % costs for deleting objects ([])
param.cutoff.cost = []; % overall cost cutoff ([] = automatic)
param.cutoff.dist = []; % cutoff for spatial distance ([] = automatic)
param.cutoff.time = -1; % cutoff for temporal distances (-1 = ignore time differences)
param.weight.dist = 1.0; % weight for distances in space (1.0)
param.weight.time.forward = 0.0; % weight for positive time distances (0.0)
param.weight.time.backward = 0.0; % weight for negative time distances (0.0)
param.weight.volume = 1000.0; % weight for distances in volume (0.0)
param.weight.intensity = 100.0; % weight for distances in intensity (0.0)
param.weight.type = Inf; % weight for different types of objects (Inf)
% printing
param.print.load = 1; % print info about loading the data files
param.print.save = 1; % print info about saving the data files
param.print.match.frames = 1; % print info about matching the frames
param.print.match.objects = 1; % print results on matching the objects
param.print.estimates = 1; % print automatic determined estimates for cutoffs
% figures
param.figure.match = 0; % generate figures for each match
param.figure.trajectories = 1; % generate figure displaying the trajectories
param.figure.stats = 1; % generate figure on statistics for the trajectories
[match, cost] = matchObjects(frame(1), frame(2), param);
figure(18)
clf
subplot(1,2,1)
plotMatchedObjects(match)
title('matches')
subplot(1,2,2)
plotMatchedCostMatrix(match, cost)
title('cost matrix')
%% full matching
[matches, traj] = trackObjects(frame, param);
%%
save('/home/ckirst/Desktop/track.mat')
