function [frames, matches, trajs] = runTracker(indir, outdir, param)
%
% [frames, matches, trajs] = runTracker(indir, outdir, param)
%
% description:
%    main routine to run the tracking by
%    - loading data in directory indir,
%    - running the tracking algrothim on all subsequent time frames, 
%    - writing analysis data to directory outdir
%
% input:
%    indir    directory of files to load object data from
%    outdir   directory of files to save data to
%    param    parameter structure containing parameter used by the various functions
%             see beginning of code for optional parameters and behavior with reduced
%             number of arguments
%             .save                save the result to folder outdir (1)
%             .min_time            start matching from this time onward ([])
%             .max_time            last possible time in trajectory of matched cells ([])
%             .max_frames          maximal number of frames to match for testing ([])
%
%
%             .figure.match        generate figures for the matching results (0)
%             .figure.trajectories generate figure for the trajectories
%             .figure.stats        generate figuer for the statistics of the trajectories
%
% output:
%    frames   array of Frame classes representing the movie
%    matches  array of Match classes containting the sequential matching information
%    trajs    array of Trajectory classes containing the trajectory information
% 
% See also: setParameter, matchFrames, findTrajectories, Frame, Match, Trajectory

if nargin < 2
   outdir = indir;
   param = [];
end

if nargin < 3
   if isstruct(outdir)
      param = outdir;
      outdir = indir;
   else
      param = [];
   end
end

sav = getParameter(param, {'save'}, 1);

min_time  = getParameter(param, {'min_time'}, []);
max_time  = getParameter(param, {'max_time'}, []);
max_frames = getParameter(param, {'max_frames'}, []);

fig_match  = getParameter(param, {'figure', 'match'}, 0);
fig_traj   = getParameter(param, {'figure', 'trajectories'}, 0);
fig_stats  = getParameter(param, {'figure', 'stats'}, 0);


%% load data
frames = loadEmbryoData(indir, param);

if ~isempty(min_time)
   times = frames.time;
   frames = frames(find(times >= min_time, 1,'first') : end);
end

if ~isempty(max_time )
   times = frames.time;
   frames = frames(1: find(times <= max_time, 1,'last'));
end

if ~isempty(max_frames)
   frames = frames(1: min(length(frames), max_frames));
end




%% match frames

if fig_match
   [matches, costs] = matchFrames(frames, param);

   for t =1:length(matches)
      figure
      clf
      subplot(1,2,1)
      plotMatchedObjects(matches(t))
      title('matches')
      subplot(1,2,2)
      plotMatchedCostMatrix(matches(t), costs{t})
      title('cost matrix')
   end
else
   matches = matchFrames(frames, param);
end   


%% find trajectories

trajs = findTrajectories(matches);

if fig_traj
   figure
   clf
   plotMatchedTrajectories(frames, trajs)
end

if fig_stats
   stats = trajs.statistics;

   figure
   subplot(1,2,1)
   hist(stats.length.values)
   title(sprintf('trajectory time lengths:\nmean:%g std:%g', stats.length.mean, stats.length.std))
   xlabel('time');

   subplot(1,2,2)
   hist(stats.dist.values)
   title(sprintf('trajectory spatial lengths:\nmean:%g std:%g', stats.dist.mean, stats.dist.std))
   xlabel('distance')
end

%% save data

if sav
   saveEmbryoData(outdir, frames, trajs, param);
end
