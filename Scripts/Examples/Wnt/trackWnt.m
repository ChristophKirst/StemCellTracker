%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Tracking Objects in 3D %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Init

close all
clear all
clc

set(0, 'DefaultFigurePosition', [1   705   560   420]);

initialize()

ijinitialize();
bfinitialize();

addpath('/home/ckirst/Science/Simulation/Matlab/StemCell3D/Scripts/Examples/Wnt')
verbose = true;

%% Load Data from a .lif file

% filename = '/home/ckirst/Science/Projects/StemCells/Experiment/Other/Wnt/wnt_clone8_again_feb09.lif';
% seriesid = 1;
% datalif = imread_bf(filename, struct('series', seriesid, 'channel', []));
% %check = imcheckimage(datalif);
% 
% 
% % last image is blank
% datalif = datalif(:,:,:,:,1:end-1);
% size(datalif)
% check = imcheckimage(datalif);
% nframes = size(datalif, ndims(datalif))
% nframes = 3;
% 
% for t = 1:nframes
%    datalif(:,:,:, 1, t) = mat2gray(datalif(:,:,:,1,t));
%    datalif(:,:,:, 2, t) = mat2gray(datalif(:,:,:,2,t));
% end
% 
% if false
%    figure(1)
%    clf
%    boxr = [1 1 1];
%    stackraw = datalif(:,:,:,1,1);
%    set(gcf, 'Name', ['Raw Stack: ' filename ' series:' num2str(seriesid) ' channel: 1']);
%    implot3d(stackraw, 'BoxRatios', boxr, 'Texture', 'z');
%    
%    
%    figure(2)
%    clf
%    stackraw = datalif(:,:,:,2,1);
%    set(gcf, 'Name', ['Raw Stack: ' filename ' series:' num2str(seriesid) ' channel: 2']);
%    implot3d(stackraw, 'BoxRatios', boxr, 'Texture', 'z');
%    
%    
%    figure(3)
%    clf
%    stackraw = datalif(:,:,:, 1, t) -0.5 * datalif(:,:,:,2,t);
%    stackraw(stackraw < 0 ) = 0;
%    set(gcf, 'Name', ['Raw Stack: ' filename ' series:' num2str(seriesid) ' channel: 1-2']);
%    implot3d(stackraw, 'BoxRatios', boxr, 'Texture', 'z');
%      
%    
%    %ijplot3d(stackraw, 'PixelDepth', 1)
% end
% 
% 
% if verbose
%    %ijplot3d(imgd, 'PixelDepth', 3);
%    t = 1;
%    lifdataf(:,:,:,1,t) = medianFilter(datalif(:,:,:,1,t),3);
%    %lifdata(:,:,:,2,t) = medianFilter(lifdata(:,:,:,2,t),3);
%    
%    lifdataf(:,:,:,1,t) = cast(256 * mat2gray(log(lifdataf(:,:,:,1,t)+eps) + 15), 'int16');
%    lifdataf(:,:,:,2,t) = cast(256 * mat2gray(log(datalif(:,:,:,2,t)+eps) + 15), 'int16');
%    
%    ijplot3d(lifdataf(:,:,:,:,1), 'PixelDepth', 5);
% 
%    t = 2;
%    lifdataf(:,:,:,1,t) = medianFilter(datalif(:,:,:,1,t),3);
% 
%    lifdataf(:,:,:,1,t) = cast(256 * mat2gray(log(lifdataf(:,:,:,1,t)+eps) + 15), 'int16');
%    lifdataf(:,:,:,2,t) = cast(256 * mat2gray(log(datalif(:,:,:,2,t)+eps) + 15), 'int16');
%    
%    ijplot3d(lifdataf(:,:,:,:,t), 'PixelDepth', 5);
%    %ijplot3d(lifdata(:,:,:,1,t), 'PixelDepth', 5);
% end


%% Load data from the Test folder

fh = FileHandler('BaseDirectory',           './Test/Data/Experiment', ...
                 'ImageDirectoryName',     '../../Images/mESCells_Wnt',...
                 'ResultDirectoryName',    '.', ...
                 'ReadImageCommandFormat', 'imread(''<file>'')',...
                 'ReadImageFileFormat',    'wnt_t<time,2>_z<z,2>.tif');

fh.info()

%%

fh.fileName('time', 1, 'z', 1)

%%
tagr = fh.fileTagRange

img = fh.readImage(struct('time', 1, 'z', 1));
[x,y] = size(img)

%% read all data at once

datalif = zeros([x, y,11, 1, 5]);
for t = tagr.time
   for z = tagr.z
      datalif(:,:,z,1,t) = fh.readImage(struct('time', t, 'z', z));
   end
end

size(datalif)


%% Segmentation

nframes = 2;
%nframes = 5;
for t = 1:nframes
   fprintf('segmenting frame %g / %g\n', t, nframes);

   %% Segment
   stack = datalif(:,:,:,1,t);
   [imgseg, stats] = segment3dWnt(stack, false);
 
   %% convert to Objects for Tracking
   objects = label2DataObjects(imgseg, datalif(:,:,:,1,t), stats, 'existing', struct('time', t, 'rescale', [1 1 5]));
   frame(t) = Frame('objects', objects);
end


%% Inspect Data Structure

imglab = frame(end).labeledImage();
figure(13)
implotlabeloutline(stack, imglab)



%% Tracking

%% Test: 2 Frames

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


%% Return Time Series

ts = TimeSeries('frames', frame, 'trajectories', traj);

