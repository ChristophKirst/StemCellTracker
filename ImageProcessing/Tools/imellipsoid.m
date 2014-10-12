function [centroids, mainaxes, varargout] = imellipsoid(labeledimage, stats)
%
% [centroids, mainaxes] = imellipsoid( labeledimage );
% [centroids, mainaxes, stats] = imellipsoid(labeledimage, stats);
%
% description:
%     finds the ellipsoids of the labeled regions in labeled images that fit the first two moments
%
% input:
%     labeledimage  labeled image 
%     stats         (optional) statistics struct of labeledimage containing PixelIdxList
%
% output
%     centroids     centroids of the ellipsoids (dim x n)
%     mainaxes      main axes of the ellipsoids (dim^2 x n)
%     stats         (optional) statistics struct of labeledimage containing additionally cacualted statistics 
%
% See also: implotellipsoid


dim = ndims(labeledimage);

if dim < 2 || dim > 3
   error('imellipsoid: expects 2d or 3d labeled image!')
end

if nargin < 2
   stats = [];
end

if dim == 2
   stats = regionprops(labeledimage, {'Centroid', 'Orientation', 'MinorAxisLength', 'MajorAxisLength'});
   
   centroids = reshape([stats.Centroid], 2,[]);
   centroids = centroids([2 1], :);
   
   mainaxes = zeros(4, length(stats));
   
   phi = - [stats.Orientation] / 180 * pi;
   
   ax = [cos(phi); sin(phi)];
   mainaxes([2 1], :) = repmat([stats.MajorAxisLength],2,1) .* ax / 2;
   
   ax = ax([2 1], :);
   mainaxes([4 3], :) = repmat([stats.MinorAxisLength],2,1) .* ax / 2;

else % dim == 3

   stats = imstatistics(labeledimage, stats,  {'Centroid', 'PixelIdxList'});
   nlabels = length(stats);
   centroids = zeros(3, nlabels);
   mainaxes = zeros(3*3, nlabels);
   isize = size(labeledimage);
   
   for i = 1:nlabels
      [cc, aa] = imellipsoid3d(isize, stats(i).PixelIdxList, stats(i).Centroid);
      centroids(:, i) = cc;
      mainaxes(:,i)  = aa;
   end

   if nargout > 2
      varargout{1} = stats;
   end
end

end



%%% helper
function [cent, mainax] = imellipsoid3d(isize, idx, cent)
      
    % extract points of the current particle
    [x, y, z] = ind2sub(isize, idx);
    
    % number of points
    n = length(idx);

    % compute approximate location of ellipsoid center
    xc = cent(1);
    yc = cent(2);
    zc = cent(3);
  
    x = (x - xc);
    y = (y - yc);
    z = (z - zc);

    points = [x y z];
        
    % principal component analysis of covariances for inertia axes
    co = cov(points) / n;
    [U, S] = svd(co);
    
    % extract length of each semi axis
    radii = 2 * sqrt(diag(S)*n)';
    
    % sort axes from greater to lower
    [radii, ind] = sort(radii, 'descend');
    
    % format U to ensure first axis points to positive x direction
    U = U(ind, :);
    if U(1,1) < 0
        U = -U;
        % keep matrix determinant positive
        U(:,3) = -U(:,3);
    end
    
    mainax = [radii(1) * U(1, :),  radii(2) * U(2,:), radii(3) * U(3,:)];
end

