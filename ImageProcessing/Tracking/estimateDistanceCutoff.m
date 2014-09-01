function cutoff = estimateDistanceCutoff(data0, data1, param)
%
% cutoff = estimateDistanceCutoff(matrix)
% cutoff = estimateDistanceCutoff(matrix, param)
% cutoff = estimateDistanceCutoff(data0, data1, param)
%
% descritpion:
%   Estimates a cutoff based on data for the distance matrix.
%   The estimate is max of the distances to the 5th nearest neighbours
%
% input:
%   matrix   the matrix of distances
%   data*    arrays of ojects with coordinate field .r
%   param    stratc as in estimateCostCutoff
%
% output:
%   cutoff   estimate for the distance cutoff
%
% See also: distanceMatrix, estimateCostCutoff

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

cutoff = estimateCostCutoff(dist, param);

end