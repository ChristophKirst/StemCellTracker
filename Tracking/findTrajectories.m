function trajetories = findTrajectories(matches)
%
% trajetories = findTrajectories(matches)
%
% description:
%   identifies trajectories as chains from the pairwise matches
%
% input:
%   matches      array of Match classes
%
% output:
%   trajetories  array of Trajectory classes
%
% See also: Match, Trajectory

% initialize

frames = {matches.objects0 matches(end).objects1 };
nframes = length(frames);


ncells = [matches.n0];
ncells(end+1) = matches(end).n1;


% find trajectories

ntraj = ncells(1);
active = 1:ntraj;
traj = num2cell(active);
times = num2cell(ones(1,ntraj));

for t = 1:nframes-1
   
   pre = cellfun(@(l) l(end), traj(active));
   post = matches(t).match(pre)';
  
   % trajectory ends
   active = active(post > 0);
   post = post(post>0);

   % trajectory continues
   traj(active) = cellfun(@(l, a) [l a], traj(active), num2cell(post),'UniformOutput', false);
   times(active) = cellfun(@(l, a) [l a], times(active), num2cell((t+1) * ones(1, length(active))),'UniformOutput', false);
   
   % trajectory starts
   new = setdiff(1:ncells(t+1), post);
   traj = [traj num2cell(new)]; %#ok<AGROW>
   times = [times num2cell((t+1) * ones(1, length(new)))]; %#ok<AGROW>
   
   active = [active, ntraj+1:length(traj)]; %#ok<AGROW>   
   ntraj = length(traj);
end
 

% construc trajectories

for n = ntraj:-1:1
   
   ti = times{n};
   tr = traj{n};
   
   for p = length(ti):-1:1
      objs(p) = frames{ti(p)}(tr(p));
   end
   
   trajetories(n) = Trajectory(objs, ti, tr);

   clear objs
end


end

