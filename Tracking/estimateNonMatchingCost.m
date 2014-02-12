function cost = estimateNonMatchingCost(data0, data1, param)
%
% cost = estimateNonMatchingCost(matrix)
% cost = estimateNonMatchingCost(matrix, param)
% cost = estimateNonMatchingCost(data0, data1, param)
%
% description:
%   Estimates the cost for creation or deletion of objects.
%   The estimate is max of the distances to the first nearest neighbours.
%
% input:
%   matrix                 the matrix of distances
%   data*                  arrays of ojects with coordinate field .r
%   param.print.estimates  print estimate
%
% output:
%   cost     estimate for the non matching cost
%
% See also: distanceMatrix


if nargin == 1
   dist = data0;
   param= [];  
elseif nargin ==2
   if isstruct(data1) || isempty(data1)
      dist = data0;
      param = data1;
   else
      dist = distanceMatrix(data0, data1);
      param = [];
   end
else
  dist = distanceMatrix(data0, data1); 
end

dist = sort(dist);
dist = dist(1,:);
cost = max(dist(dist < Inf));


if getParameter(param, {'print', 'estimates'}, 0)
   fprintf('estimateNonMatchingCost: estimated non/linking cost: %g\n', cost);
end

end