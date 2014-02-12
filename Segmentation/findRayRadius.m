function [radius, rays, r0, polx, poly] = findRayRadius(image, imagegrad, seeds, param)
%
% [radius, rays, r0, polx, poly] = findRayRadius(image, imagegrad, seeds, param)
%
% description: 
%    starts rays from r0 and stops them if certain criteria are matched
%    returns position on rays and the corresponding polygon coordinates
%
% input:
%    image      instensity image to be segmented
%    imagegrad  gradinet of the intensity image
%    seeds      bw image of the size of image with white pixels indicating the sseds for the rays
%    param      struct with parameter entries
%               .nrays                              number of rays (15)
%               .cutoff.radius                      maximal radius (30)
%               .center_averaging                   radius for dermining average intensity and gradient at seed centers (0 = center pixel)
%               .threshold.background               intensities below this value are backgournd (0 = none)
%               .threshold.relative_change          stop propagation if reltative change is larger (0.5, Inf = none)
%               .threshold.absolute_change          stop propagation if absolute change is larger (Inf = none)
%               .threshold.relative_gradient_peak   minimal peak height in gradient profile (Inf = none)
%               .threshold.absolute_gradient_peak   minimal peak height in gradient profile (Inf = none)
%               .trough.drop                        minimal trough falloff in itensity profile (Inf = do not use trough method)
%               .trough.max_rise                    stop trough finding if profile rises more than this
%               .trough.rise                        minimal subsequent trough increase in itensity profile
%               .trough.noise                       if no increase reset trough minima by ignoring fluctuations less than this              
%               .smooth                             use robust LOWESS smoothing on the radius outline (false)
%               .plot                               generate figures with detailed results (false)
%
% output:
%    radius     distances along the rays of the form [r1_1, r1_2, ...; r2_1, r2_2, ...; ...]
%    rays       single set of ray directions along the seeds are extened
%    r0         starting points of the rays
%    polx,poly  x and y positions of the estimated points of the estimated polygons [x1_1, x1_2, ...; x2_1, x2_2, ...;...]
%
% note: 
%    all coordinates in pixel coordinates h,w
%
% See also: segmentByRays, segmentByActiveRays

% parameter 
 
nrays         = getParameter(param, {'nrays'}, 15);

cutoff_radius = getParameter(param, {'cutoff', 'radius'}, 30);

averaging_radius = getParameter(param, {'center_averaging'}, 0);  % radius for dermining average intensity and gradient at seed centers

background             = getParameter(param, {'threshold', 'background', },   0.0);  % definately stop ray below this intensity
relative_change        = getParameter(param, {'threshold', 'relative_change'}, 0.8); % threshold for relative change in intensity
absolute_change        = getParameter(param, {'threshold', 'absolute_change'}, Inf); % threshold for relative change in intensity
absolute_gradient_peak = getParameter(param, {'threshold', 'absolute_gradient_peak'},   Inf); % minimal peak height in gradient profile
relative_gradient_peak = getParameter(param, {'threshold', 'relative_gradient_peak'}, Inf); % minimal peak height in gradient profile

trough_max_rise   = getParameter(param, {'trough', 'max_rise'},   0.1);    % stop trough finding if profile rises more than this
trough_drop       = getParameter(param, {'trough', 'drop'},   Inf);        % minimal trough falloff in itensity profile
trough_rise       = getParameter(param, {'trough', 'rise'},   0.1);        % minimal subsequent trough increase in itensity profile
trough_noise      = getParameter(param, {'trough', 'noise'},   0.1);       % if no increase reset trough minima by ignoring fluctuations less than this

smooth_radius     = getParameter(param, {'smooth'},   false);             % smooth final radii with spline

show              = getParameter(param, {'plot'},   false);             % plot final results


over_sample  = 2;
nraypixel = over_sample * cutoff_radius; % oversampling by factor of 2 makes it smoother

% initialize

[r0x, r0y] = find(seeds);
r0 = [r0x'; r0y'];

phi = 0:2*pi/nrays: (2*pi - 2*pi/nrays);
rays = [cos(phi); sin(phi)];

nr0s = size(r0x,1);

radius = cutoff_radius * ones(nr0s, nrays);
if nargout > 3
   polx = zeros(nr0s, narays);
   poly = zeros(nr0s, narays);
end


% loop over seeds

for ir0 = 1:nr0s

   % center rays at r0
   r0i = repmat(r0(:, ir0), 1, nrays);   
   r1 = r0i + cutoff_radius * rays;
   y = [r0i(1,:)', r1(1,:)']; % r0 etc are given in pixel coordiantes
   x = [r0i(2,:)', r1(2,:)']; % [h,w] = [y,x]

   % intensity profiles 
   profile = zeros(nraypixel, nrays);
   for i = 1:nrays
      profile(:, i) = improfile(image, x(i,:), y(i,:), nraypixel); %note: improfile assumes spatial coordinates !!!
   end

   %
   % find estimate for positions on poly
   %

   % bound by backgournd
   if background > 0 
      for i=nrays:-1:1
         rm = find(profile(:,i) < background, 1, 'first');
         if ~isempty(rm)
            radius_background(i) = rm;
         else
            radius_background(i) = nraypixel;
         end
      end
   end

   
   if absolute_change < Inf || relative_change < Inf || intensity_trough < Inf    
      if averaging_radius == 0
         center_intensity = profile(1,:);
      else
         center_intensity = double(imfiltervalues(image, r0i(:,1), averaging_radius));
         center_intensity = repmat(center_intensity, 1, nrays);       
      end
   end
   
      
   %absolute change
   if absolute_change < Inf
      
      absolute_change_profile = abs(profile - repmat(center_intensity, nraypixel,1) );

      for i = nrays:-1:1
         ac = find(absolute_change_profile(:,i) > absolute_change, 1, 'first');
         if ~isempty(ac)
            radius_absolute_change(i) = ac;
         else
            radius_absolute_change(i) = nraypixel;
         end
      end
   end

   % relative change in intensity profile
   if relative_change < Inf || trough_drop < Inf    
      relative_profile = profile ./ repmat(center_intensity, nraypixel,1);
      relative_change_profile = abs(relative_profile - 1);
   end
   
   if relative_change < Inf   
      for i = nrays:-1:1
         rc = find(relative_change_profile(:,i) > relative_change, 1, 'first');
         if ~isempty(rc)
            radius_relative_change(i) = rc;
         else
            radius_relative_change(i) = nraypixel;
         end
      end
   end

   
   %trough or fall off in relative profile larger than trough_drop and 
   if trough_drop < Inf      
      for i = nrays:-1:1
         radius_trough(i) = findTroughReset(relative_profile(:,i), trough_max_rise, trough_drop, trough_rise, trough_noise);
      end
   end
    
   
   %peak in gradient
   if absolute_gradient_peak < Inf || relative_gradient_peak < Inf   
      gradient_profile = zeros(nraypixel, nrays);
      for i = 1:nrays
         gradient_profile(:, i) = improfile(imagegrad, x(i,:), y(i,:), nraypixel);
      end
   end
      
   if absolute_gradient_peak < Inf 
      for i = nrays:-1:1
         gp = findPeaks(gradient_profile(:,i), absolute_gradient_peak, 1);
         if ~isempty(gp)
            radius_absolute_gradient_peak(i) = gp(1);
         else
            radius_absolute_gradient_peak(i) = nraypixel;
         end
      end
   end
   
   
   % peak in relative gradient profile
   if relative_gradient_peak < Inf
      if averaging_radius == 0
         center_gradient = gradient_profile(1,:);
      else
         center_gradient = imfiltervalues(imagegrad, seeds, averaging_radius);
      end

      relative_gradient_profile = gradient_profile ./ repmat(center_gradient, nraypixel,1);
      
      % check peaks in gradient
      for i = nrays:-1:1
         gp = findPeaks(relative_gradient_profile(:,i), relative_gradient_peak, 1);
         if ~isempty(gp)
            radius_relative_gradient_peak(i) = gp(1);
         else
            radius_relative_gradient_peak(i) = nraypixel;
         end
      end 
   end

   
   %
   % choose radius
   %
   
   radius_min = nraypixel * ones(1, nrays);
   if background > 0 
      indx = radius_min > radius_background;
      radius_min(indx) = radius_background(indx);
   end
   if absolute_change < Inf
      indx = radius_min > radius_absolute_change;
      radius_min(indx) = radius_absolute_change(indx);
   end  
   if relative_change < Inf
      indx = radius_min > radius_relative_change;
      radius_min(indx) = radius_relative_change(indx);
   end  
   if trough_drop < Inf
      indx = radius_min > radius_trough;
      radius_min(indx) = radius_trough(indx);
   end
   if absolute_gradient_peak < Inf
      indx = radius_min > radius_absolute_gradient_peak;
      radius_min(indx) = radius_absolute_gradient_peak(indx);
   end  
   if relative_gradient_peak < Inf
      indx = radius_min > radius_relative_gradient_peak;
      radius_min(indx) = radius_relative_gradient_peak(indx);
   end  
   
   radius_min(radius_min > nraypixel) = nraypixel;
   
   if smooth_radius
      radius_smooth = smoothLowess(radius_min, 3);
      radius(ir0, :) = round(radius_smooth / over_sample); % correct for over_sampling
   else
      radius(ir0, :) = round(radius_min / over_sample); % correct for over_sampling
   end

   if nargout > 3
      pos = vertcat(radius(ir0, :),radius(ir0, :)) .* rays + r0;
      polx(ir0, :) = pos(1,:);
      poly(ir0, :) = pos(2,:);
   end
   
   %
   % debug / plotting
   %
   
   if show
      figure
      plotRays(image, imagegrad, rays, r0(:, ir0), radius(ir0, :), radius_min, cutoff_radius)

      if background == 0; 
         radius_background = [];
      end 
      if absolute_change == Inf && relative_change == Inf && intensity_trough == Inf
         center_intensity = [];
      end
      if absolute_change == Inf
         absolute_change_profile = [];
         radius_absolute_change = [];
      end  
      if relative_change == Inf && trough_drop == Inf    
         relative_profile = [];
         relative_change_profile = [];
      end 
      if relative_change == Inf
         radius_relative_change = [];
      end  
      if trough_drop == Inf
         radius_trough = [];
      end
      if  absolute_gradient_peak == Inf && relative_gradient_peak == Inf
         gradient_profile = [];
      end
      if absolute_gradient_peak == Inf
         radius_absolute_gradient_peak = [];
      end  
      if relative_gradient_peak == Inf
         center_gradient = [];
         relative_gradient_profile = [];
         radius_relative_gradient_peak = []; 
      end  

      figure
      plotProfiles(profile, absolute_change_profile, relative_profile, relative_change_profile, gradient_profile, relative_gradient_profile,...
                   center_intensity, center_gradient,...
                   radius_min, radius_background, radius_absolute_change, radius_relative_change, radius_trough, radius_absolute_gradient_peak, radius_relative_gradient_peak);
   end
             
end % loop ove seeds

end






%% Subroutines 

%
% trough with reset
%
% finds a dropoff of at least drop, if no such dropoff returns length(x) 
% if x rises above rise again it returns position of the minimum
% if x does not rise again resets to first point in x in which change is larger then noise level (< drop!)
% we allow for an initial peak not larger than rise
function trough = findTroughReset(x, rise1, drop, rise2, reset)

   n = numel(x);
   i = 0; d = 0; xu = x(1); xd = x(1);

   while i < n

      i = i + 1; xi = x(i);

      switch d    
         case -1
            if xi < xd % still decaying
               trough = i; xd = xi;
            elseif xi >= xd + rise2 % strong enough growth -> found trough
               return
            end
         case 0 % for initial state
            if xi <= xu - drop % trend is decaying
               d = -1; trough = i;
            elseif xi >= xd + rise1; % rise1 threshold reached -> return
               trough = i;
               return
            end
            if xi < xd  % update max / min
               xd = xi;
            elseif xi > xu
               xu = xi;
            end
      end
   end

   
   % did not succeed, check if we need to rest after drop
   if d == -1
      xt = x(trough);
      r = find(abs(x(trough-1:-1:1) - xt) > reset, 1, 'first');
      if ~isempty(r)
         trough = trough - r + 2;
      end
   else
      trough = n;
   end

end




%% visualization for debuggin

% show rays and ploygons on images
function plotRays(image, imagegrad, rays, r0, radius, radius_min, cutoff_radius)

   % draw rays on images 
   nrays = length(radius);
   pos =  cutoff_radius * rays + repmat(r0, 1, nrays);   
   y = [repmat(r0(1), nrays, 1), pos(1,:)'];  % we use pixel coordinates
   x = [repmat(r0(2), nrays, 1), pos(2,:)'];  % [h,w] = [y,x]
   pos =  [radius; radius] .* rays + repmat(r0, 1, nrays); 
   posm = [radius_min; radius_min] .* rays + repmat(r0, 1, nrays); 
   
   col = linspace(0,1, nrays+1);
   col = col(1:end-1);
   col = hsv2rgb([col', ones( nrays,1), ones(nrays,1)]);

   ax(1)= imsubplot(1,2,1);
   imshow(image);
   for i=1:nrays
      line(x(i,:), y(i,:), 'Color', col(i,:))
   end
   linecl(posm(2,:), posm(1,:), 'Color', 'g');
   linecl(pos(2,:), pos(1,:), 'Color', 'r');

   
   ax(2) = imsubplot(1,2,2);
   imshow(imagegrad);
   for i=1:nrays
      line(x(i,:), y(i,:), 'Color', col(i,:))
   end
   linecl(posm(2,:), posm(1,:), 'Color', 'g');
   linecl(pos(2,:), pos(1,:), 'Color', 'r');
      
   linkaxes(ax, 'xy')
end
   
   
% plot profiles along rays
function plotProfiles(profile, absolute_change_profile, relative_profile, relative_change_profile, gradient_profile, relative_gradient_profile,...
                center_intentisty, center_gradient,...
                radius_min, radius_background, radius_absolute_change, radius_relative_change, radius_trough, radius_absolute_gradient_peak, radius_relative_gradient_peak) 
  
nrays = size(profile, 2);
sp = size(profile);

col = linspace(0,1, nrays+1);
col = col(1:end-1);
col = hsv2rgb([col', ones( nrays,1), ones(nrays,1)]);
set(gcf,'DefaultAxesColorOrder',col)

ax(1) = subplot(4,2,1);
axi = 1;
plot(profile)
hold on
if ~isempty(radius_min)
   plot(radius_min, profile(sub2ind(sp, radius_min,1:nrays)), '*');
end
if ~isempty(radius_background)
   plot(radius_background, profile(sub2ind(sp, radius_background,1:nrays)), 'b+');
end
if ~isempty(center_intentisty)
   plot(ones(1,nrays),center_intentisty, 'gx')
end

title('intensity profile,  *min, +background, x center intensity')

if ~isempty(relative_profile)
   axi = axi +1;
   ax(axi) = subplot(4,2,3);
   plot(relative_profile)
   hold on
   if ~isempty(radius_trough)
      plot(radius_trough, relative_profile(sub2ind(sp, radius_trough,1:nrays)), '*');
   end
   title('relative intensity,  *trough')
end

if ~isempty(relative_change_profile)
   axi = axi +1;
   ax(axi) = subplot(4,2,5);
   plot(relative_change_profile)
   hold on
   if ~isempty(radius_relative_change)
      plot(radius_relative_change, relative_change_profile(sub2ind(sp, radius_relative_change,1:nrays)), '*');
   end
   title('relative intensity change,  *rel change')
end

if ~isempty(absolute_change_profile)
   axi = axi +1;
   ax(axi) = subplot(4,2,7);
   plot(absolute_change_profile)
   hold on
   if ~isempty(radius_absolute_change)
      plot(radius_absolute_change, absolute_change_profile(sub2ind(sp, radius_absolute_change,1:nrays)), '*');
   end
   title('absolute intensity change,  *abs change')
end

if ~isempty(gradient_profile)
   axi = axi +1;
   ax(axi) = subplot(4,2,2);
   plot(gradient_profile)
   hold on
   if ~isempty(radius_absolute_gradient_peak)
      plot(radius_absolute_gradient_peak, gradient_profile(sub2ind(sp, radius_absolute_gradient_peak,1:nrays)), '*');
   end
   if ~isempty(center_gradient)
      plot(ones(1,nrays),center_gradient, 'b+') 
   end
   title('gradient, *abs grad peak, +center grad')
end

if ~isempty(relative_gradient_profile)
   axi = axi +1;
   ax(axi) = subplot(4,2,4);
   plot(relative_gradient_profile)
   hold on
   if ~isempty(radius_relative_gradient_peak)
      plot(radius_relative_gradient_peak, relative_gradient_profile(sub2ind(sp, radius_relative_gradient_peak,1:nrays)), '*');
   end
   title('relative gradient, *abs grad peak')
end

linkaxes(ax, 'x');

end


