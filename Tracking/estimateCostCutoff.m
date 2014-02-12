function cutoff = estimateCostCutoff(cost, param)
%
% cutoff = estimateCostCutoff(cost, param) 
%
% descritpion:
%   Estimates a cutoff based on the cost matrix.
%   The estimate is max of the cost to the 5th nearest neighbours
%
% input:
%   cost                    the matrix of linking costs
%   param.print.estimates   switch to print the estimate
%
% output:
%   cutoff  estimate for the cost cutoff
%
% See also: costMatrixObjectMatching

if nargin < 2
   param = [];
end

n = min(5, size(cost,1));

cost(cost == Inf) = 0;
cost = sort(cost);
cutoff = max(cost(n,:));

if getParameter(param, {'print', 'estimates'}, 0)
   fprintf('estimateCostCutoff: estimated cost cutoff: %g\n', cutoff);
end

end