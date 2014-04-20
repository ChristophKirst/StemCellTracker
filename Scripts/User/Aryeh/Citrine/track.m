function track()
%
% skript for tracking the Citrine data
%

%%
if false
   %%
   apath = '/home/ckirst/Science/Simulation/Matlab/StemCell3D/Scripts/User/Aryeh/Citrine';
   addpath(apath)
   
   close all
   clear all
   clc

   set(0, 'DefaultFigurePosition', [1   705   560   420]);

   verbose = true;
   initialize()
end



%% Data Format
[fileformat, range] = tagformat('/home/ckirst/Media/ChristophsData/Science/Projects/StemCells/Experiment/Other/Aryeh/Citrine/p1', {'z', 't'});
fileformat

zrange = unique(range(1,:));
nzslices = length(zrange);
fprintf('zrange: %s\n', var2char([min(zrange), max(zrange)]))

trange = unique(range(2,:));
fprintf('trange: %s\n', var2char([min(trange), max(trange)]))

nframes = length(trange);

%% Ceate Frames
clear frame
nframes = 10;
tstart = 4;
for t = tstart:(tstart + nframes-1)
   %% Load Image

   fileformat = '/home/ckirst/Media/ChristophsData/Science/Projects/StemCells/Experiment/Other/Aryeh/Citrine/p0/LSM_31_2014_04_13__11_05_51_z<z,2>_t<t,3>_p0.tif';
   filename = tags2name(fileformat, [3, 1]);
 
   xrange = 500:600;
   yrange = 900:1000;

   %t = 1;
   img = zeros(length(xrange), length(yrange), length(zrange));
   for z = 1:nzslices
      filename = tags2name(fileformat, [zrange(z), t]);
      img(:,:,z) = imread(filename, 'PixelRegion',{ [xrange(1), xrange(end)], [yrange(1), yrange(end)]});
   end
   
   figure(1)
   colormap(gray)
   implottiling(img)
   
   %%
   img = double(img)/2^16;  
   size(img)   
   
   
   %% Segmentation
   
   [imgseg, stats] = segment3d(img, true, 2 * t);
   
   %%
   %size(imgseg)
   %size(img)
   %max(imgseg(:))
   %length(unique(imgseg(:)))
   nsegs = length(stats)

   %imgr = imrelabel(imgseg);
   %max(imgr(:))
   %length(unique(imgr(:)))
   

   %% Create Objects for Tracking
   [objects, stats] = label2Objects(imgseg, img, stats, setParameter('time', t));
   
   frame(t-tstart+1) = Frame('objects', objects);
end



%% test track

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
param.optimize = true;       % optimize match using an additional optimal coordinate transformation and second matching

% cost matrix
param.cost.creation  = [];   % costs for creating new objects ([])
param.cost.deletion = [];    % costs for deleting objects ([])

param.cutoff.cost = [];      % overall cost cutoff ([] = automatic)
param.cutoff.dist = [];      % cutoff for spatial distance ([] = automatic)
param.cutoff.time = -1;      % cutoff for temporal distances (-1 = ignore time differences)

param.weight.dist = 1.0;          % weight for distances in space (1.0)

param.weight.time.forward = 0.0;  % weight for positive time distances (0.0)
param.weight.time.backward = 0.0; % weight for negative time distances (0.0)
param.weight.volume = 1000.0;     % weight for distances in volume (0.0)
param.weight.intensity = 100.0;   % weight for distances in intensity (0.0)
param.weight.type = Inf;          % weight for different types of objects (Inf)

% printing
param.print.load = 1;           % print info about loading the data files
param.print.save = 1;           % print info about saving the data files

param.print.match.frames  = 1;  % print info about matching the frames
param.print.match.objects = 1;  % print results on matching the objects

param.print.estimates = 1;      % print automatic determined estimates for cutoffs

% figures
param.figure.match = 0;          % generate figures for each match
param.figure.trajectories = 1;   % generate figure displaying the trajectories
param.figure.stats = 1;          % generate figure on statistics for the trajectories


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









