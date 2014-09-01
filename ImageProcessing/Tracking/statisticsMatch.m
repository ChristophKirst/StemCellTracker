function stats = statisticsMatch(match, cost)
%
% stat = statisticsMatch(match)
% stat = statisticsMatch(match, cost)
%
% description:
%    returns statistics on the match given the cost matrix cost
%
% input:
%    match   match class or array of matches
%    data*   arrays of objects
%    cost    cost matrix
%
% See also: Match

if ~isa(match, 'Match')
   error('expect match to be of class Match')
end
   
if nargin < 2
   cost = match.cost;
end

% spatial distances
[r0, r1] = match.toCoordinates();
stats.dist.values = sum((r0- r1).^2, 1);
stats.dist.mean = mean(stats.dist.values);
stats.dist.std = std(stats.dist.values,1);


% matching cost 
if ~isempty(cost)
   pairs = match.toArray(-1,-1);
   pairs = pairs(min(pairs,[],2)>0, :);
   indx = sub2ind(size(cost), pairs(:,1), pairs(:,2));
   stats.cost.values = cost(indx);
   stats.cost.mean = mean(stats.cost.values);
   stats.cost.std = std(stats.cost.values,1);
end
