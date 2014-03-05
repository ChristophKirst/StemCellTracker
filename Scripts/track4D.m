%% Tracking in 4D

% Init

close all
clear all
clc

set(0, 'DefaultFigurePosition', [1   705   560   420]);

initialize()

ijinitialize();

%% Data files


datadir = '/var/run/media/ckirst/Seagate Backup Plus Drive/Data/Movie/';

nframes = length(dir(fullfile(datadir, 'W1F188T*Z01C1.tif')))

nframes = 5;


for f = 1:nframes
   filenames{f} = fullfile(datadir, ['W1F188T' num2str0(f, 4) '*C1.tif']);
end

%%


%%

nfiles = length(filenames);
nframes = nfiles; 
frame0 = 1;

for f = nframes:-1:frame0
   [seg, prop] = segment4D(filenames{f});
   
   cent = {prop.Centroid};
   cent = cellfun(@(x) x([2,1,3]), cent, 'UniformOutput', false);
   [prop.Centroid] = cent{:};

   segdata{f,2} = prop;
   segdata{f,1} = seg;
   
end

%% save the stuff

save('./Save/TrackingTest.mat', 'segdata');

%%

load('./Save/TrackingTest.mat');

nframes = length(segdata); 
frame0 = 1;


%% convert segmentation to Frames / Objects

for f = frame0:nframes
   segprop = segdata{f,2};
   no = length(segprop);

   clear objects;
   for i = no:-1:1
      objects(i) = Object(i, f, segprop(i).Centroid', segprop(i).Area, median(segprop(i).PixelValues));
   end
   
   frame(f) = Frame(filenames{f}, objects);
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
param.weight.volume = 0.0;        % weight for distances in volume (0.0)
param.weight.intensity = 0.0;     % weight for distances in intensity (0.0)
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


[matches, traj] = trackObjects(frame, param);





