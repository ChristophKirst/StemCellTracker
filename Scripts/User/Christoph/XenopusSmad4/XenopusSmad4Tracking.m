%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Segment and Track Xenopus Smad 4 movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clear classes
close all
clc
initialize
bfinitialize

addpath('./Scripts/User/Christoph/XenopusSmad4');

verbose = true;

initializeParallelProcessing(12) % number of processors



%% Initialize Files

fns = '/data/Science/Projects/StemCells/Experiment/Xenopus/Smad4/12-17-14 venus smad4 mutant movie2/wild2.lsm';
is = ImageSourceBF(fns);
%is.setCellDataFormat('XYZ', 'T');
clc
is.printInfo

% C=1 nucleus, c=2 smad4

%% 3D Segmentation

for time = 1:is.dataSize('T')
   imgraw = is.data('T', time, 'C', 1);
   
   %% max projection
   imgmp = max(imgraw, [], 3);
   
   
   %% find nuclei in max proj
   [imgseg, imgmask, imgmed] = XenopusSmad4Segment2D(imgmp, verbose, time);

   
   %% crate objects
   [objects, stats] = label2Objects(imgseg, imgmp, 'time', time);
   
   frame(time) = Frame('objects', objects, 't', time);
   
   % save image data
   imgs{time} = {imgmp, imgseg};
   
end


%% Tracking 

%method
param.optimize = false;       % optimize match using an additional optimal coordinate transformation and second matching

% cost matrix
param.cost.creation  = [];   % costs for creating new objects ([])
param.cost.deletion = [];    % costs for deleting objects ([])

param.cutoff.cost = [];      % overall cost cutoff ([] = automatic)
param.cutoff.dist = [];      % cutoff for spatial distance ([] = automatic)
param.cutoff.time = -1;      % cutoff for temporal distances (-1 = ignore time differences)

param.weight.dist = 1.0;          % weight for distances in space (1.0)

param.weight.time.forward = 0.0;  % weight for positive time distances (0.0)
param.weight.time.backward = 0.0; % weight for negative time distances (0.0)
param.weight.volume = 100000.0;     % weight for distances in volume (0.0)
param.weight.intensity = 1000.0;   % weight for distances in intensity (0.0)
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



%%

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


%% Relabel Full Trajectories

% only keep full length traj
ids = {traj.objids};
ll = cellfun(@length, ids);
pos = ll == is.dataSize('T');

trajfull = traj(pos);


ntraj = length(trajfull);

imgtraj = cell(maxtime,1);
statstraj = cell(maxtime,1);
for t = 1:maxtime
   imgtraj{t} = imrelabel(imgs{t}{2});
   statstraj{t} = imstatistics(imgtraj{t}, 'PixelIdxList');
   imgtraj{t} = zeros(size(imgs{t}{2}));
end

for l = 1:ntraj
   tids = trajfull(l).frameids;
   objs = trajfull(l).objids;
 
   for i = 1:length(tids);
      t = tids(i);
      o = objs(i);
      
      imgt = imgtraj{t};
      statst = statstraj{t};
      imgt(statst(o).PixelIdxList) = l;
      
      imgtraj{t} = imgt;
   end
end


for t = 1:maxtime
   imgt = imgtraj{t};
   imgt(1,1) = ntraj + 10;
   imgtraj{t} = imgt;
end

%% Save Results for ImageJ

fexp = '/home/ckirst/Desktop/XenopusSmad4_WT/Label_T<T,3>.tif';

maxtime = l   imgraw = is.data('T', time, 'C', 1);
   
   %% max projection
   imgmp = max(imgraw, [], 3);
   ength(imgs);
cols = colorcube(ntraj+10);

movie = zeros(is.dataSize('X'), is.dataSize('Y'), 1, 3, is.dataSize('T'));

for t = 1:maxtime
   imgovl = imoverlaylabel(imgs{t}{1}, impixelsurface(imgtraj{t}), false, 'color.map', cols,  'color.shuffle', 'noshuffle');
   %implottiling(imgovl, 'tiling', [4,3]);

   imwrite(imgovl, tagExpressionToString(fexp, 'T', t));
      
   if verbose
         figure(9); clf
         implot(imgovl);
         drawnow;
   end

   movie(:,:,1,:,t) = imgovl;
   
end

%%
ijinitialize
ijplot5d(movie)



%%

% segmentation in 2d + t ??

img2dt = zeros(is.dataSize('X'), is.dataSize('Y'), is.dataSize('T'));
for t = 1:12
   
   imgraw = is.data('T', t, 'C', 1);
   imgmp = max(imgraw, [], 3);

   img2dt(:,:,t) = imgmp;
end

%%
ijplot3d(img2dt, 'PixelDepth', 10)





