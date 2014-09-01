function plotMatchedCostMatrix(match, cost)
%
% plotMatchedCostMatrix(match)
% plotMatchedCostMatrix(match, cost)
%
% deescription:
%    plots the cost matrix and indicates the matches 
%
% input:
%    match    Match object
%    cost     cost matrix if Match does not contain one
%
% See also: Match

if ~isa(match, 'Match')
   error('expect match to be of class Match')
end

if nargin < 2
   cost = match.cost;
end

pairs = match.toArray([], []);

hold on
imagesc(cost);
scatter(pairs(:,2),pairs(:,1), 50, [1 1 1]);
axis([0 match.n1 0 match.n0] + 0.5);

xlabel('post'); ylabel('pre');

hold off
