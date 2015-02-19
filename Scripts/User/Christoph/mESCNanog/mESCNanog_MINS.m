%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Mouse Nanog visualize MINs results %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clear classes
close all
clc
initialize
bfinitialize
ijinitialize

verbose = true;

addpath('./Scripts/User/Christoph/mESCNanog');
addpath('./Interface/User/EmbryoData');


initializeParallelProcessing(12) % number of processors

tiling = [5,4];

%% Data

datadir = '/data/Science/Projects/StemCells/Experiment/Mouse/Nanog/';
%dexp = 'T<F,3>/T<F,3>.tif';
%fulldexp = fullfile(datadir, dexp);
fulldexp = fullfile(datadir, '12Aug14FGFonNanogH2BGFP-700_movie.lsm');
dns = tagExpressionToFiles(fulldexp);

isd = ImageSourceBF(fulldexp);
isd.printInfo

%% MINs Results
% note: z stack here is saved as time T !
resultdir = '/data/Science/Projects/StemCells/Experiment/Mouse/Nanog/12Aug14FGFonNanogH2BGFP-700_MINS_1/';
fexp = '12Aug14FGFonNanogH2BGFP-700_movie_channel=0001_frame=<F,4>_segmentation.tiff';

fullfexp = fullfile(resultdir , fexp);
fns = tagExpressionToFiles(fullfexp)

is = ImageSourceFiles(fullfexp);


%% Plot Data
f = 1;

img = isd.data('T', f, 'C', 1);
imgseg = is.data('F', f);

figure(1); clf
implottiling(imoverlaylabel(img, impixelsurface(imgseg), false), 'tiling', tiling)


%% Plot 5D data

nT = 3;
imgmovie = zeros([is.dataSize, 3, nT]);

for t = 1:nT
   
   img = isd.data('T', t, 'C', 1);
   imgseg = is.data('F', t);
   imgmovie(:,:,:,:, t) = imoverlaylabel(img, impixelsurface(imgseg), false);
end

ijplot5d(imgmovie)



%% Track MINS result
clc

param.load.min = 2;       % at least one object in frame
param.load.change = 0.2;  % at most 20% change in number of objects

param.print.load = true;
param.print.match.objects = true;
param.print.match.optimization = false;

param.optimize = true;

%vns = dir([resultdir, '*.csv']);
%frm = loadEmbryoDataFile([resultdir, vns(1).name])
 
%%
 
frames = loadEmbryoData(resultdir, param);

for t = 1:length(frames)
   
   dat = frames(t).r;
   %indx = dat(2,:) > 340;
   %frames(t).objects = frames(t).objects(indx);
   
end

[matches, costs] = matchFrames(frames, param);


%% plot the result 
% 
% for t =1:length(matches)
%    figure
%    clf
%    subplot(1,2,1)
%    plotMatchedObjects(matches(t))
%    title('matches')
%    subplot(1,2,2)
%    plotMatchedCostMatrix(matches(t), costs{t})
%    title('cost matrix')
% end


%% Determine trajectories

traj = findTrajectories(matches);

figure(17)
clf
plotMatchedTrajectories(frames, traj)


%% Plot some statistics
lt = [traj.duration];
trajL = traj(lt == 67);

figure(18)
clf
plotMatchedTrajectories(frames, trajL)


%% 

figure(19); clf; hold on
for k = 1:length(trajL)
   i = [trajL(k).objects.intensity];
   plot(i)
end


%% Plot in 5d

moviedir = fullfile(resultdir, 'tracking');
mkdir(moviedir)

outfile = fullfile(moviedir, 'T<T,3>Z<Z,2>.tiff');

ids = {traj.objids};
ll = cellfun(@length, ids);

% only keep full length traj
%pos = ll == is.dataSize('F');

tmax = is.cellSize('F');
ntraj = length(traj);

cols = colorcube(ntraj);
cols = cols(randperm(size(cols,1)),:);

imgmovie = zeros([is.dataSize, 3, tmax]);

for t = 1:tmax
   
   % load data
   img = isd.data('T', t, 'C', 1);
   imgseg = is.data('F', t);
   
   % determine color map
   [tpos, topos, tobjs] = traj.frameSlice(t);
   
   %maps is [tobjs.id] -> tpos
   colm = zeros(length(tpos), 3);
   colm([tobjs.id], :) = cols(tpos, :);

   imgmovie(:,:,:,:, t) = imoverlaylabel(mat2gray(img), impixelsurface(imgseg), false, 'color.map', colm,  'color.shuffle', 'noshuffle');
   
   imgwrt = imfrmtReformat(imgmovie(:,:,:,:,t), 'XYZC', 'XYCZ');
   
   for z = 1:size(imgwrt,4)
      imwrite(imgwrt(:,:,:,z), tagExpressionToString(outfile, 'T', t, 'Z', z));
   end
end


%%
ijplot5d(moviedir)



















%%
%% some statistics

stats = statisticsTrajectory(trajs);

figure
subplot(1,2,1)
hist(stats.length.values)
title(sprintf('trajectory time lengths:\nmean:%g std:%g', stats.length.mean, stats.length.std))
xlabel('time');

subplot(1,2,2)
hist(stats.dist.values)
title(sprintf('trajectory spatial lengths:\nmean:%g std:%g', stats.dist.mean, stats.dist.std))
xlabel('distance')


%% saving data

saveEmbryoData('./Test/Out', frames, trajs)



%% run full Tracker 

runTracker('./Test/Data', './Test/Out')


%% run full Tacker with standard parameter and a filter 

path(path, './Test')

param = setParameter();

param.filter =  @testFilter;

runTracker('./Test/Data', './Test/Out', param)


%% run full Tacker with test parameter 
path(path, './Test')

param = setParameterTest();

param.filter = @testFilter;

[frames, matches, trajs] = runTracker('./Test/Data.all', './Test/Out', param);
















