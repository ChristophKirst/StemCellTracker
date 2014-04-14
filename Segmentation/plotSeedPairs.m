function plotSeedPairs(imglab, pairs, col)
%
% plotSeedPairs(imglab, pairs)
%
% description: plot lines between pairs
%
% input:
%    imglab   labeled image with single seed for each index
%    pairs    pairs of indices pair(i,1) <-> pair(i,2)

if nargin < 3
   col = 'b';
end

np = length(pairs);

pos = imstatistics(imglab, 'Centroid');
pos = round([pos.Centroid]);

hold on
for i = 1:np  
   xy =[pos(:, pairs(i,1)), pos(:, pairs(i,2))];
   line(xy(1,:), xy(2,:), 'Color', col)
end

end