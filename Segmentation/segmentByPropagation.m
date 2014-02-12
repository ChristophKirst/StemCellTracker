function [segments, distances] = segmentByPropagation(image, label, mask, lambda, radius_dist, radius_center)
%
% [segments, distances] = segmentByPropagation(label, image, mask, lambda)
%
% description: 
%     segments image using propagation algorithm to 8 nearest neighbours
%     neighbours are included if distance is small enough
%     distance is calucalted as relative intensity change to make algorith
%     robust against differently illuminated objects
%
% input:
%    image         instensity image to be segmented
%    lables        starting lables
%    mask          restrict propagation to this mask
%    lambda        weight for spatial distances (wieght for intensity distances = 1)
%    radius_dist   radius of square to calculate intensity differences 
%    radius_center radius of square to calcualte intensity of seeds for normalization (-1 for no normalization) 
%
% output:
%    segments   labels for segments
%    distances  distances to initial seeds
%


if nargin < 3
   mask = ones(size(image));
end
if nargin < 4
   lambda = 1.0;
end
if nargin < 5
   radius_dist = 3;
end
if nargin < 6
   radius_center = 3;
end

[segments, distances] = segmentByPropagationMEX(image, label, mask, lambda, radius_dist, radius_center);

end

   