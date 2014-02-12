function stats = statisticsTrajectory(traj)
%
% stat = statisticsTrajectory(traj)
%
% description:
%    returns some useful statistics on the trajectories traj
%
% See also: Trajectory

if ~isa(traj, 'Trajectory')
   error('statisticsTrajectory: expect traj to be of class Trajectory')
end

% discrete trajectory lengths
stats.n = length(traj);
stats.length.values = cellfun(@length, {traj.objects});
stats.length.mean = mean(stats.length.values);
stats.length.std = std(stats.length.values,1);


% spatial distances

ntraj = length(traj);
stats.dist.values(ntraj) = 0;

for i=1:ntraj 
   xyzt = traj(i).r;
   stats.dist.values(i) = sum(sqrt(sum((xyzt(:, 2:end) - xyzt(:, 1:end-1)).^ 2)));
end   

stats.dist.mean = mean(stats.dist.values);
stats.dist.std = std(stats.dist.values,1);


