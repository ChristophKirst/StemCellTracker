function [stacksegmented, segmentprops] = segment4dWntfromLeicaCK(filename, seriesid)
% Segmentation in 4D

%% Init

close all
clear all
clc

set(0, 'DefaultFigurePosition', [1   705   560   420]);

initialize()

ijinitialize();

addpath('/home/ckirst/Science/Simulation/Matlab/StemCell3D/Scripts/User/Eric/Wnt')
verbose = true;

%% Load Data

filename = '/home/ckirst/Science/Simulation/Matlab/StemCell3D/Test/Images/Develop/Wnt/wnt_clone8_again_feb09.lif';
seriesid = 1;
datalif = imread_bf(filename, struct('series', seriesid, 'channel', []));
%check = imcheckimage(datalif);


% last image is blank
datalif = datalif(:,:,:,:,1:end-1);
size(datalif)
check = imcheckimage(datalif);
nframes = size(datalif, ndims(datalif))

for t = 1:nframes
   datalif(:,:,:, 1, t) = mat2gray(datalif(:,:,:,1,t));
   datalif(:,:,:, 2, t) = mat2gray(datalif(:,:,:,2,t));
end

if false
   figure(1)
   clf
   boxr = [1 1 1];
   stackraw = datalif(:,:,:,1,1);
   set(gcf, 'Name', ['Raw Stack: ' filename ' series:' num2str(seriesid) ' channel: 1']);
   implot3d(stackraw, 'BoxRatios', boxr, 'Texture', 'z');
   
   
   figure(2)
   clf
   stackraw = datalif(:,:,:,2,1);
   set(gcf, 'Name', ['Raw Stack: ' filename ' series:' num2str(seriesid) ' channel: 2']);
   implot3d(stackraw, 'BoxRatios', boxr, 'Texture', 'z');
   
   
   figure(3)
   clf
   stackraw = datalif(:,:,:, 1, t) -0.5 * datalif(:,:,:,2,t);
   stackraw(stackraw < 0 ) = 0;
   set(gcf, 'Name', ['Raw Stack: ' filename ' series:' num2str(seriesid) ' channel: 1-2']);
   implot3d(stackraw, 'BoxRatios', boxr, 'Texture', 'z');
     
   
   %ijplot3d(stackraw, 'PixelDepth', 1)
end


if verbose
      %ijplot3d(imgd, 'PixelDepth', 3);
   t = 1;
   lifdataf(:,:,:,1,t) = medianFilter(datalif(:,:,:,1,t),3);
   %lifdata(:,:,:,2,t) = medianFilter(lifdata(:,:,:,2,t),3);
   
   lifdataf(:,:,:,1,t) = cast(256 * mat2gray(log(lifdataf(:,:,:,1,t)+eps) + 15), 'int16');
   lifdataf(:,:,:,2,t) = cast(256 * mat2gray(log(datalif(:,:,:,2,t)+eps) + 15), 'int16');
   
   ijplot3d(lifdataf(:,:,:,:,1), 'PixelDepth', 5);
   
   
   t = 2;
   lifdataf(:,:,:,1,t) = medianFilter(datalif(:,:,:,1,t),3);

   lifdataf(:,:,:,1,t) = cast(256 * mat2gray(log(lifdataf(:,:,:,1,t)+eps) + 15), 'int16');
   lifdataf(:,:,:,2,t) = cast(256 * mat2gray(log(datalif(:,:,:,2,t)+eps) + 15), 'int16');
   
   ijplot3d(lifdataf(:,:,:,:,t), 'PixelDepth', 5);
   
   
   %ijplot3d(lifdata(:,:,:,1,t), 'PixelDepth', 5);
end



%% Segmentation

dataseg = zeros(size(datalif));
for t = 1:nframes
   %stack = datalif(:,:,:, 1, t) - 0.5 * datalif(:,:,:,2,t);
   %stack(stack < 0) = 0;
   stack = datalif(:,:,:,1,t);
   dataseg(:,:,:,1,t) = segment3dWnt(stack, false);
end


%% Convert to Objects and add statistics

for t = 1:nframes
   objects = label2objects(datalif(:,:,:,1,t), dataseg(:,:,:,1,t), t, [1 1 5]);
   frame(t) = Frame(filename, objects);
end


%% Tracking



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









