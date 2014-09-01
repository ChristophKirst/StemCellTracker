function [match, cost] = matchFrames(frames, param)
%
% [match, cost] = matchFrames(frames, param)
% 
% description:
%    Compute the optimum matching between subsequent frames in data
%
% input: 
%    frames   array of Frame objects
%    param    parameter struct with entries as in 
%             costMatrixObjectMatching, costMatrixFrameMatching
%             for negative parameter .cutoff.dist, .cost.deletion, .cost.creation
%             estimates for the cutoffs are calculated only once    
%             .print.match.frames  print info on the result of the frame matching
%                          
% output:
%    match    array of Match objects
%    cost     cell array of cost matrices
%    
% See also:  costMatrixObjectMatching, matchFrames, Match, Frame

% initialize

if nargin < 2
   param = [];
end

print_match = getParameter(param, {'print', 'match', 'frames'}, 1);



% for negative parameter we calculate the estimates once
% use [] for automatic detection in each round with little extra computation

cutoff_dist = getParameter(param, {'cutoff', 'dist'}, []);
cost_creation = getParameter(param, {'cost', 'creation'}, []);
cost_deletion = getParameter(param, {'cost', 'creation'}, []);

if ((~isempty(cutoff_dist) && cutoff_dist < 0) ...
   ||(~isempty(cost_creation) && cost_creation < 0) ...
   ||(~isempty(cost_deletion) &&  cost_deletion < 0)) && lenght(data) > 1

   dist = distanceMatrix(data(1), data(2));

   if (cutoff_dist < 0)
      param.cutoff.dist = estimateDistanceCutoff(dist);
   end
   if (cost_creation < 0)
      param.cost.creation = estimateNonLinkingCost(dist);
   end
   if (cost_deletion < 0)
      param.cost.deletion = estimateNonLinkingCost(dist);
   end
end


% find best linking between successive frames

nframes = length(frames);
match(nframes-1) = Match();
if (nargout > 1)
   cost{nframes-1} = [];
end

for t = 1:nframes-1
  
   if print_match
      fprintf('matchFrames: matching step %d / %d, matching frames %d -> %d with times %g -> %g\n', t, nframes-1,t,t+1, frames(t).time, frames(t+1).time);
   end
   
   [match(t), co] = matchObjects(frames(t), frames(t+1),  param);

   if (nargout > 1)
      cost{t} = co;
   end
end


end
  