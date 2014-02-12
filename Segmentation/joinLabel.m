function reduced = joinLabel(image, gradient, label, param)
%
% reduced = joinLabel(image, label, param)
%
% description:
%    If intensity changes and gradients along the conencting line of two
%    labels is indicating the same cell then join label.
%
% input:
%    image    original image
%    gradient gradient of the original image
%    label    label of objects
%    param    parameter struct with entries
%             .threshold.background    % below this intensity objects end definately
%             .threshold.change        % maximal rel change in intensitiy above objects are assumed to be different
%             .threshold.gradient      % maximal absolute gradient change above objects are joint 
%             .cutoff.distancen        % maximal distance between labels (= 20)
%             .reference_radius        % radius of center object to calculate reference mean intensity (=3)
%
% output:
%    rdeuced reduced labels
%
% See also: findPeaks


% initialize

if nargin < 4
   param = [];
end

cutoff_distance    = getParameter(param, {'cutoff', 'distance'}, 20);
threshold_change   = getParameter(param, {'threshold', 'change'}, 0.5);
threshold_gradient = getParameter(param, {'threshold', 'gradient'}, 0.1);
threshold_backgournd =  getParameter(param, {'threshold', 'background'}, 0.0);
reference_radius =  getParameter(param, {'refernece_radius'}, 3);


% find all pairs within distance cutoff

[x, y] = find(label > 0);
xy = [x, y]'
dist = distanceMatrix(xy);

figure
imshow(mat2gray(dist))


dist(tril(ones(size(dist)))) = Inf;
[ind0, ind1] = find(dist < cutoff_distance);

% generate intensity and gradient profiles along the pairs

npairs = length(ind0);
dist = dist(sub2ind(size(dist), ind0, ind1));

xy0 = xy(:, ind0);
xy1 = xy(:, ind1);
x01 = [xy0(1, :); xy1(1,:)];
y01 = [xy0(2, :); xy1(2,:)];

xy


means = imfiltervalues(image, xy, reference_radius);
means0 = means(ind0);
means1 = means(ind1);

% visualize this for debugging

figure 
imshow(imoverlay(image, label > 0))
for p = 1:nparis
   line(x01(:,p), y01(:,p))
end

return



for p = npairs:-1:1
   profile{p} = improfile(image, x01(:, p), y01(:, p), round(dist(p))+1);
   profile{p} = profile{p} / min(means0(p), means1(p));
 
  % gradprofile{p} = improfile(image, x01(:, p), y01(:, p), roudn(dist(p)+1);
end



% calcualte statistics



% join labels




end

