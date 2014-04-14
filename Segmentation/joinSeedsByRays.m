function [imgjoin, varargout] = joinSeedsByRays(imglab, img, imggrad, param)
%
% imgjoin = joinSeedsByRays(img, imglab, param)
% imgjoin = joinSeedsByRays(img, imggrad, imglab, param)
% [imgjoin, pairs, joins] = joinSeedsByRays(...)
%
% description:
%    If intensity changes and gradients along the connecting line of two
%    labels do not change sufficiently then the seeds are joined and
%    the conencting line is added to the labeled pixel
%
% input:
%    img      intensity image
%    imggrad  gradient of img
%    imglabel labeled or bw image
%    param    parameter struct with entries
%             .threshold.min           if profile falls below this intensity objects are different (0)
%             .threshold.max           if profile stais above thus absolute intensity objects are joined (inf)
%             .threshold.change        maximal rel change in intensitiy above objects are assumed to be different 
%             .threshold.gradient      maximal absolute gradient change below objects are joined 
%             .cutoff.distance         maximal distance between labels (= 20)
%             .averaging.ksize         ksize to calculate reference mean intensity (=3)
%             .addline                 add a line between joined label (true)
%
% output:
%    imgjoin  reduced joined seeds image
%    pairs    (optional) all considered pairs as rows
%    joins    (optional) joined pairs
% 
% todo: need to extend improfile for 3d images !
%
% See also: findPeaks

% initialize

if nargin < 3
   param = [];
elseif nargin < 4
   if isstruct(imggrad)
      param = imggrad;
      imggrad = [];
   else
      param = [];
   end
end
 
cutoff_distance      = getParameter(param, {'cutoff', 'distance'}, 15);
threshold_change     = getParameter(param, {'threshold', 'change'}, 1.5);
threshold_gradient   = getParameter(param, {'threshold', 'gradient'}, 0.1);
threshold_min        = getParameter(param, {'threshold', 'min'}, 0.0);
threshold_max        = getParameter(param, {'threshold', 'max'}, inf);
averaging_ksize      = getParameter(param, {'averaging', 'ksize'}, 3);
addline              = getParameter(param, {'addline'}, 1);

% find center points
centr = imstatistics(imglab, {'Centroid', 'PixelIdxList'});
pixidx = {centr.PixelIdxList};
centr = round([centr.Centroid]);

isize = size(imglab);
for d = 1:2
   ch = centr(:,d);
   ch(ch > isize(d)) = isize(d);
   centr(:,d) = ch;
end
centr(centr <= 0) = 1;

%clable = zeros(isize);
%clabel(imsub2ind(isize, cent')) = 1;
%[pairs, dist] = imlabellocalpairs(clabel, cutoff_distance, 'distance');

% find all pairs within distance cutoff
dist = distanceMatrix(centr);
[i, j] = find(dist < cutoff_distance & dist > 0);   
pairs = sort([i, j], 2);
pairs = unique(pairs, 'rows');
dist = dist(sub2ind(size(dist), pairs(:,1), pairs(:,2)));
npairs = length(pairs);

if nargout > 1
   varargout{1} = pairs;
end


%calcualte local averages for relevant points
pairsidx = unique(pairs);
means0 = zeros(1, max(pairsidx));
means0(pairsidx) = double(imfiltervalues(img, centr(:, pairsidx), averaging_ksize));

% find pairs to join
join = [];
checkgrad = ~isempty(imggrad);
jp = false;

for p = 1:npairs
   p1 = pairs(p,1); p2 = pairs(p,2);
   c1 = centr(:,p1); c2 = centr(:, p2);
   x = [c1(1), c2(1)]; y = [c1(2), c2(2)];
   
   profile =  improfile(img, y, x, 2*round(dist(p)));
   
   if p < 10
   figure(p + 30)
   subplot(1,2,1)
   plot(profile)
   end
   
   % lover absolute threshold -> if we cross background dont join
   if min(profile) < threshold_min
      continue
   end

   % high intensity joins -> if profile is above threshold definately join
   if threshold_max < inf
      if min(profile) > threshold_max
         jp = true;
      end
   end
   
   if ~jp
      % chang in rel intensity
      profile = profile / min(means0(p1), means0(p2));
   
      if p < 10
      subplot(1,2,2)
      plot(profile)
      end
      
      if max(profile)-min(profile) > threshold_change
         continue
      end
   
      % check gradient profile
      if checkgrad
         gradprofile = improfile(imggrad, y, x, 2*round(dist(p)));
         if max(gradprofile) < threshold_gradient
            continue
         end
      end
   else
      jp = false;
   end

   % all checks succesfull -> join piars 
   join = [join; pairs(p,:)]; %#ok<AGROW>
         
end

%join

% finally join labels
imgjoin = imglab;

if nargout > 2
   varargout{2} = join;
end


nj = size(join, 1);
if nj == 0
   return
end

pid1 = join(:,1);

for p = 1:nj
   
   %relabel the exsiting label
   imgjoin(pixidx{join(p,2)}) = pid1(p);
   pid1(pid1 == join(p,2)) = pid1(p);
  
   % add pixel line
   if addline
      imgjoin = impixelline(imgjoin, centr(:, join(p,1)), centr(:, join(p,2)), pid1(p));
      %size(imgjoin)
   end 
end


end




