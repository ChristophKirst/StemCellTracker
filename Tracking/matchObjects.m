function [match, cost] = matchObjects(data0, data1, param)
%
% [match, cost] = matchObjects(data0, data1, param)
%
% description:
%   finds an optimal matching between two arrays of objects data*
%
% input:
%   data*   classes with field objects or array of objects to be matched
%           fields for the objects in of array: r, time, volume, intensity, identity
%
%   param   is parameter structure with entries as in costMatrixObjectMatching
%           .optimize                   optimizes the matching by finding a suitable coordinate transformation ([])
%           .print.match.optimization   print the optimization result
%           .print.match.objects        print info on the result of the matching
%
% output:
%   match   Match class
%   cost    used cost matrix for matching
%
% See also: Match, costMatrixObjectMatching, optimalAssociationMatrix, optimalTransformation

% initialize
if nargin < 3
   param = [];
end

optimize     = getParameter(param, {'optimize'}, []);
print_match  = getParameter(param, {'print', 'match', 'objects'}, 0);
print_opt    = getParameter(param, {'print', 'match', 'optimization'}, 0);

if isentry(data0, 'objects')
   data0 = data0.objects;
end
if isentry(data1, 'objects')
   data1 = data1.objects;
end


% match 
[match, cost] = findMatch(data0, data1, param);

if ~isempty(optimize) && optimize > 0
   % optimal coord transformation
   [X0, X1] = match.toCoordinates();

   [R, T, C] = optimalTransformation(X0,X1);
   
   if print_opt
      disp 'matchObjects: optimal transformation:'
      R %#ok<NOPRT>
      T %#ok<NOPRT>
      C %#ok<NOPRT>
   end
   
   data0t = Object(data0); % deep copy of relevant Object data only!
   data0t = data0t.transformCoordinates(R, T, C);

   % match transformed data
   [match, cost] = findMatch(data0t, data1, param);
   
   % restore original objects
   match.objects0 = data0;
end


if print_match
   
   fprintf('matchObjects: trying to match %d to %d objects, %d matches made', match.n0, match.n1, sum(match.match>0))
   if ~isempty(optimize) && optimize > 0
      fprintf(' using optimization\n\n');
   else
      fprintf('\n\n');
   end
      
end

end




function [match, cost] = findMatch(data0, data1, param)

% find best linking
cost = costMatrixObjectMatching(data0, data1, param);
A = optimalAssociationMatrix(cost);

% convert matrix to permutation list

osize = size(cost,1) - 1;
nsize = size(cost,2) - 1;

[iLink, jLink] = find(A(1:osize,:));
jLink(jLink==(nsize+1)) = -1;
[~, perm] = sort(iLink);
jLink = jLink(perm);

match = Match(jLink, data0, data1);

end


