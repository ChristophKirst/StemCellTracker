%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Mouse Nanog visualize MINs results                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialization 

%% Initialize

clear all
clear classes
close all
clc
initialize
bfinitialize
ijinitialize
addpath('./Interface/Other/MINs');

verbose = true;

% put path of this script here
addpath('./Scripts/User/Christoph/mESCNanog');

initializeParallelProcessing(12) % number of processors

verbose = true;  % plot intermediate results

%% Setup Data Sources
clc

basedir = '/data/Science/Projects/StemCells/Experiment/Mouse/Nanog/12Aug14FGFonNanogH2BGFP-700';

% Microscope Data
lsmfile = fullfile(basedir, '12Aug14FGFonNanogH2BGFP-700_movie.lsm');
isd = ImageSourceBF(lsmfile);
%isd.setCaching(cache_data);
isd.printInfo

% MINs Data
% note: z stack here is saved as time T !
minsdir = fullfile(basedir, '/12Aug14FGFonNanogH2BGFP-700_MINS_1/');
minsfileexps = '12Aug14FGFonNanogH2BGFP-700_movie_channel=0001_frame=<F,4>_segmentation.tiff';
minsfileexps = fullfile(minsdir , minsfileexps);
minsfiles    = tagExpressionToFiles(minsfileexps);

isMINS = ImageSourceFiles(minsfileexps);
%isMINS.setCaching(cache_mins);
nT = isMINS.cellSize;
nZ = isMINS.dataSize('T');
nX = isMINS.dataSize('X');
nY = isMINS.dataSize('Y');

fprintf('found %d segmented images in folder %s\n', nT, minsdir);

%% Plot Data

if verbose
   f = 1;
   img = isd.data('T', f, 'C', 1);
   imgseg = isMINS.data('F', f);
   
   figure(1); clf
   implottiling(imoverlaylabel(img, impixelsurface(imgseg), false), 'tiling', [5,4])
end

%% Plot 5D data

if verbose 
   nt = 2;
   imgall = zeros([isMINS.dataSize, 3, nt]);
   
   for t = 1:nt
      img = isd.data('T', t, 'C', 1);
      imgseg = isMINS.data('F', t);
      imgall(:,:,:,:, t) = imoverlaylabel(img, impixelsurface(imgseg), false);
   end
   
   ijplot5d(imgall, 'write', true)
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Tracking

%% Setup Tracking Parameter
clc

% tracking parameter
param.load.min = 2;       % at least one object in frame
param.load.change = 0.2;  % at most 20% change in number of objects to regard next frame as valid

% method
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
param.weight.volume = 0.5;        % weight for distances in volume (0.0)
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
param.figure.trajectories = 0;   % generate figure displaying the trajectories
param.figure.stats = 0;          % generate figure on statistics for the trajectories

 
%% Track MINS result

% load
frames = loadEmbryoData(minsdir, param);

% potential filter on frames 
% for t = 1:length(frames)
%    dat = frames(t).r;
%    %indx = dat(2,:) > 340;
%    %frames(t).objects = frames(t).objects(indx);
% end

% match
[matches, costs] = matchFrames(frames, param);


%% Plot Track Result 

if verbose 
   %tend = length(matches);
   tend = 1;
   for t =1:tend
      figure
      clf
      subplot(1,2,1)
      plotMatchedObjects(matches(t))
      title('matches')
      subplot(1,2,2)
      plotMatchedCostMatrix(matches(t), costs{t})
      title('cost matrix')
   end
end


%% Determine Trajectories

traj = findTrajectories(matches);
fprintf('found %d trajectories!\n', length(traj));

%%
if verbose
   figure(17); clf
   plotMatchedTrajectories(frames, traj)
end


%% Basic Trajectory Statistics

stats = statisticsTrajectory(traj);

figure(10)
subplot(1,2,1)
hist(stats.length.values)
title(sprintf('trajectory time lengths:\nmean:%g std:%g', stats.length.mean, stats.length.std))
xlabel('time');

subplot(1,2,2)
hist(stats.dist.values)
title(sprintf('trajectory spatial lengths:\nmean:%g std:%g', stats.dist.mean, stats.dist.std))
xlabel('distance')


%% Plot Trajectories of Certain Duration

trajDuration = [traj.duration];
trajFull = traj(trajDuration == nT);
fprintf('%d / %d trajecotries are full length\n', length(trajFull), length(traj))


%%
figure(18); clf
plotMatchedTrajectories(frames, trajFull)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Analysis / Visualization

%% Setup colors etc

ntraj = length(traj);

cols = colorcube(ntraj+1); cols = cols(2:end, :); 
cols = cols(randperm(size(cols,1)),:);
size(cols)

%% Load Data 

% set ranges
trange = 1:3;
zrange = 1:nZ;
xrange = 1:nX;
yrange = 1:nY;

nts = length(trange); nzs = length(zrange);

% load the data
imgdata = isd.data('C', 1, 'T', trange, 'Z', zrange, 'X', xrange, 'Y', yrange);
imgdata = mat2gray(imgdata);
isize = size(imgdata);
ndataxyz = prod(isize(1:3));

% read label via region props
imglabel = cell(1, nts);
for ti = 1:nts
   t = trange(ti);
   img = isMINS.data('F', t);
   img = img(xrange, yrange, zrange);
   img = impixelsurface(img);
   
   % trajectory label
   imgstats = regionprops(img, 'PixelIdxList');
   
   % convert label to full 5d indices
   pixids = {imgstats.PixelIdxList};
   
   offset = ndataxyz * (ti -1) * 3;
   offsetR = offset;
   offsetG = offset + ndataxyz;
   offsetB = offset + 2*ndataxyz;
   
   for i = 1:length(pixids)
      pixids{i} = {pixids{i} + offsetR, pixids{i} + offsetG, pixids{i} + offsetB};
   end
   
   imglabel{t} = pixids;
end

%group label according to trajectories
trajlabel = cell(1,ntraj);
for i = 1:ntraj
  
   fids = traj(i).frameids;
   ids = ismember(fids, trange);
   fids  = fids(ids);

   oids = traj(i).objids;
   oids = oids(ids);
   
   labR = []; labG = []; labB = [];
   for ti = 1:length(fids);
      t = fids(ti);
      oid = oids(ti);
      labids = imglabel{t}{oid};
      labR = [labR; labids{1}]; labG = [labG; labids{2}]; labB = [labB; labids{3}]; %#ok<AGROW>
   end
   trajlabel{i} = {labR, labG, labB};
end


%% Plot 5d data 


% Select Trajectories
tids = find([traj.duration] == nT);

imgplot = repmat(imgdata, 1, 1, 1, 3, 1);
%for ti = tids
for ti = tids
   imgplot(trajlabel{ti}{1})  = cols(ti, 1);
   imgplot(trajlabel{ti}{2})  = cols(ti, 2);
   imgplot(trajlabel{ti}{3})  = cols(ti, 3);
end

ijplot5d(imgplot, 'write', true)

%% Plot Statistics 

figure(24); clf;
subplot(1,2,1); hold on
for k = tids
   i = [traj(k).objects.intensity];
   plot(i, 'Color', cols(k, :))
end

subplot(1,2,2); hold on
for k = tids
   i = [traj(k).objects.volume];
   plot(i, 'Color', cols(k, :))
end



%% Exmaple: Trajectorys starting with high expression  Select Trajectories

[tri, ~, objs] = traj.frameSlice(1);
tids = intersect(tri([objs.intensity] > 20), find([traj.duration] == nT))


%% Plot 5d data 

imgplot = repmat(imgdata, 1, 1, 1, 3, 1);
for ti = tids
   imgplot(trajlabel{ti}{1})  = cols(ti, 1);
   imgplot(trajlabel{ti}{2})  = cols(ti, 2);
   imgplot(trajlabel{ti}{3})  = cols(ti, 3);
end

ijplot5d(imgplot, 'write', true)


%% Plot Statistics 

figure(24); clf;
subplot(1,2,1); hold on
for k = tids
   i = [traj(k).objects.intensity];
   plot(i, 'Color', cols(k, :))
end

subplot(1,2,2); hold on
for k = tids
   i = [traj(k).objects.volume];
   plot(i, 'Color', cols(k, :))
end


%% Saving trajectory data for excel

% saveEmbryoData('basedir/result.csv', frames, trajs)








