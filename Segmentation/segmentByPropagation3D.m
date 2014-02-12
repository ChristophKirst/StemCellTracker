function [segments, distances] = segmentByPropagation3D(image, label, mask, param)
%
% [segments, distances] = segmentByPropagation3D(label, image, mask, param)
%
% description: 
%     segments image using propagation algorithm to 26 nearest neighbours
%     neighbours are included if distance is small enough
%     distance is calucalted as relative intensity change to make algorith
%     robust against differently illuminated objects
%
% input:
%    image         instensity image to be segmented
%    lables        starting lables / pre segmented / local max
%    mask          restrict propagation to this mask
%    param         parameter struct
%                  .weight.distance     weight for spatial distanace
%                  .weight.intensity    weight for relatie intensity difference
%                  .weight.area         weight for area (0 = ignore)
%                  .weight.surface      weight for surface (0 = ignore)
%                  .weight.convexity    weight for penalizing non convex shapes
%                  .target.area         traget area
%                  .target.surface      target surface
%                  .background          intensities below this value are considered as background
%                  .hsize.average       size of average filter for calculating the relative intensity differences
%                  .hsize.center        size of average filter to calculate reference center intensity
%
% output:
%    segments   labels for segments
%    distances  distances to initial seeds
%

if nargin < 3
   mask = ones(size(image));
end
if nargin < 4
   param = [];
end

par.weight_distance = getParameter(param, {'weight', 'distance'}, 1);
par.weight_intensity = getParameter(param, {'weight', 'intensity'}, 1);

par.weight_radius = getParameter(param, {'weight', 'radius'}, 0);
par.target_radius = getParameter(param, {'target', 'radius'}, 5);
par.target_radius2 = par.target_radius^2;

par.weight_surface = getParameter(param, {'weight', 'surface'}, 0);
par.target_surface = getParameter(param, {'target', 'surface'}, round(4 * pi * target_radius^2));

par.weight_volume = getParameter(param, {'weight', 'volume'}, 0);
par.target_volume = getParameter(param, {'target', 'volume'}, round(4/3 * pi * target_radius^3));

par.weight_convexity = getParameter(param, {'weight', 'convexity'}, 10);

par.backgournd = getParameter(param, {'background'}, 0);

par.hsize_av = getParameter(param, {'hsize', 'average'}, 10);
if length(par.hsize_av) < 3
   par.hsize_av = [par.hsize_av par.hsize_av par.hsize_av];
else
   par.hsize_av = par.hsize_av(1:3);
end
[par.offsets_lo, par.offsets_hi] = filteroffsets(par.hsize_av);
par.size_av = prod(par.hsize_av);
if par.size_av <= 0
   error('segmentByPropagation3D: averaging of intensity needs to be over at least one pixel!')
end

% initialize objects 

pos = find(label>0);



pixel_list = h







end


% distance between pixel and segmented object
function d = distance(image, pixel, obj, param)
  
   d = 0;
   r = pixel.r;
   
   % space
   if param.weight_space > 0 || param.weight_radius > 0
      dist2 = sum((pixel.r - obj.r).^2);
      d = d + param.weight_space * dist2;
     
      %if dist2 > param.target_radius2 
         d = d + param.weight_radius * (sqrt(dist2) - param.target_radius);
      %end
   end
   
   % intensity
   if param.weight_intensity > 0
      
      int = image((param.offsets_lo(1) + r(1)) : (param.offsets_hi(1) + r(1)),...
                  (param.offsets_lo(2) + r(2)) : (param.offsets_lo(2) + r(2)),...
                  (param.offsets_lo(3) + r(3)) : (param.offsets_lo(3) + r(3)));
      int = sum(int(:)) / param.size_av;
      
      d = d + param.weight_intensity * (int / obj.intensity - 1)^2;
   end

   % area
   %if param.weight_area > 0
   %   d = d + weight_area  * (obj.area + 1 - param.target_area);
   %end
      
   % volume
   if param.weight_volume > 0
      d = d + weight_volume  * (obj.volume + 1 - param.target_volume);
   end

   % convexity
   if param.weight_convexity > 0
      
      diff = abs(obj.hull - obj.mask);
      d = d + param.weight_convexity * sum(diff(:)); %one should only incoorporate the change due to the pixel !!!
     
   end


end


function obj = add_to_object(obj, pixel)
   r = pixel.r;
   obj.mask(r(1), r(2), r(3)) = 1;
   %obj.volume = sum(obj.mask(:));  
   obj.volume = obj.volume + 1;
   
   %obj.area
   obj.hull = bwconvhull3d(obj.mask);
   
end


