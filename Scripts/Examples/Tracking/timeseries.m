%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Tracking diluted cells in time %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clear classes
close all
clc

addpath('/home/ckirst/Science/Simulation/Matlab/StemCell3D/Scripts/Examples/Tracking')

initialize

%% File Handler
exp = Experiment('name', 'Test', 'description', 'Test some TimeTracking',...
                 'BaseDirectory',          './Test/Data/Experiment', ...
                 'ImageDirectoryName',     '../../Images/hESCells_Cytoo',...
                 'ResultDirectoryName',    '.', ...
                 'ReadImageCommandFormat', 'imread(''<file>'')',...
                 'ReadImageFileFormat',    'Big1_CFP_<time>.tif');

exp.info()



%%
nframes = 5;
verbose = false;

for t = 3:nframes
   %%
   img = fh.readImage('time', t);
   [imgseg, stats] = segment2dDilute(img, verbose, 0);
   
   %% Create Objects and Frame form labeled image

   param = setParameter('time' ,  0, ...   % time for objects (0)
                        'rescale',1, ...   % rescale coordinates r by this factor ([1, 1(, 1)])
                        'method', 'median'); % how to calcualte the intensity in Object, a string of any function, 'none' = dont calcualte ('median')

   objs = label2Objects(imgseg, img, stats, param);
   
   frame(t) = Frame('objects', objs, 't', t);
   
   clear stats
end

%%
frame(1:2)= [];


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
param.weight.volume = 1000.0;        % weight for distances in volume (0.0)
param.weight.intensity = 100.0;     % weight for distances in intensity (0.0)
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


%% Plot object moving in time

imgtraj = traj(1).labeledImage();
figure(42); clf
implot(imgtraj)


%% Generate Time Series

ts = TimeSeries('frames', frame, 'trajectories', traj);


%% Save 
exp.result = ts;
exp.saveExperiment('tracking.mat');



%% plot some statistics

figure(42); clf; 
subplot(1,2,1)
hist(double([frame(1).objects.intensity]), 15)
subplot(1,2,2)
hist(double([frame(1).objects.volume]), 15)

