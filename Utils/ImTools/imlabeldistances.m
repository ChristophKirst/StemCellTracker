function dist = imlabeldistances(label)
%
% pairs = imlabeldistances(label)
%
% description:
%      determines the pairwise Euclidian dinstances between the labeled points
%
% input:
%      label     labeled image
% 
% output:
%      dist      distance matrix

idx = find(label);

r0 = imind2sub(size(label), idx)';

dist = distanceMatrix(r0);

end
