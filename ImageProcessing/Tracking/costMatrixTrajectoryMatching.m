function cost = costMatrixTrajectoryMatching(data, traj, param)
%
% cost = costMatrixTrajectoryMatching(data, traj, param)
%
% description:
%   calcalates the cost matrix for linking / gap closing / mergin / slippting 
%   trajectories
%
% input:
%   data      array of TrackingTimeFrame data objects
%   traj
%   param     parameter structure with entries:
%             .cost.start           - cost for starting a trajectory
%             .cost.end             - cost for starting a trajectory
%             .cost.merge           - cost for creating a merge
%             .cost.split           - cost for creating a split
%             .cost.angle           - cost for deviations in angle
%            
%             .cutoff.cost          - overall cost cutoff
%             .cutoff.dist          - cutoff for spatial distance
%             .cutoff.time_forward  - cutoff for time to look for possible links in the future
%             .cutoff.time_backward - cutoff for time to look for possible links in the past 
%
% output:
%   cost      the cost matrix, last row/colum represent creation / deletion cost
%
% todo:
%             distance to the image borders / relinking bad matches /
%             motion propagation
%
% reference:
%             
%
% See also: distanceMatrix, costMatrixObjectMatching


error('not implemented yet')

% initialize

cost_start = param.cost.start;
cost_end   = param.cost.end;
cost_merge = param.cost.merge;
cost_split = param.cost.split;
cost_angle = param.cost.angle;

cutoff_cost= param.cutoff.cost;
cutoff_dist =param.cutoff.dist;
cutoff_time_forward = param.cutoff.time_forward;
cutoff_time_backward = param.cutoff.time_forward;

nframes = length(data);


%% linking cost between trajectory start and end points

startIds = traj.startIds();
endIds = traj.endIds();

% data points

dim = data.dim;

for i = startIds
   



%% mergin / splitting cost












% distance cost

if nargin < 5 || isempty(dist_cutoff) 
   % automatic detection
   cost(1:n, 1:m) = distanceMatrix(data0, data1);
   dist_cutoff = estimateDistanceCutoff(cost(1:n, 1:m));
   cost(cost > dist_cutoff) = Inf;
elseif dist_cutoff < 0
   % no cut_off
   cost(1:n, 1:m) = distanceMatrix(data0, data1);
else
   % given distance cutoff
   cost(1:n, 1:m) = distanceMatrix(data0, data1, dist_cutoff);
end


% object creation and deletion

if nargin < 3 || isempty(creation_cost)
   creation_cost = estimateNonLinkingCost(cost(1:n, 1:m));
end
if nargin < 4 || isempty(deletion_cost)
   deletion_cost = creation_cost;
end

for i = 1:n
   cost(i,m+1) = deletion_cost;
end
for j = 1:m
   cost(n+1,j) = creation_cost;
end
cost(n+1,m+1) = Inf;

end  