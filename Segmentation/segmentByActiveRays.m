function segments = segmentByActiveRays(image, imagegrad, seeds, mask, param)
%
% segments = segmentByActiveRays(image, param)
%
% description: 
%     segments image by sending active rays from seeds
%
% input:
%    image      instensity image to be segmented
%    seeds      bw image containnig the seeds
%    mask       restrict segmentation to this mask
%    param      struct with parameter entries
%
% output:
%    segments   labeled image indicating the segmented objects
%

%%%
%%% initialize 
%%%

%
% parameter
%

% rays / geometry
nrays = getParameter(param, {'nrays'}, 15);
cutoff_radius = getParameter(param,{'cutoff', 'radius'}, 20); % max radius

% energy evolution
delta_radius = getParameter(param, {'delta'}, 1.0);    % magnitude of change of radii in single time step
time_scale = getParameter(param, {'time_scale'}, 1.0); % overall time scale for evolution of radii
max_steps  = getParameter(param, {'max_steps'}, 100); % maximal number of time steps
stop_precision  = getParameter(param, {'stop_precision'}, 1.0); % maximal number of time steps

% energy functional 
par.weight_angle      = getParameter(param, {'energy', 'weight', 'angle'}, 1.0);
par.weight_surface    = getParameter(param, {'energy', 'weight', 'surface'}, 0.0);
par.weight_area       = getParameter(param, {'energy', 'weight', 'area'}, 1.0);
par.weight_radius     = getParameter(param, {'energy', 'weight', 'radius'}, 0.0);
par.weight_max_radius = getParameter(param, {'energy', 'weight', 'max_radius'}, 100.0);
par.weight_backgroudn = getParameter(param, {'energy', 'weight', 'background'}, 100.0);

par.weight_gradient =  getParameter(param, {'energy', 'weight', 'gradient'}, 10.0);
par.weight_intensity = getParameter(param, {'energy', 'weight', 'intensity'}, 0.0);
par.weight_variation = getParameter(param, {'energy', 'weight', 'variation'}, 10.0);

par.target_radius  = getParameter(param, {'energy', 'target', 'radius'}, 10.0);
par.target_surface = getParameter(param, {'energy', 'target', 'surface'}, 10.0);
par.target_area    = getParameter(param, {'energy', 'target', 'area'}, 10.0);
par.target_angle   = getParameter(param, {'energy', 'target', 'angle'}, (1-2/nrays) * pi); % perfect n-gon inner angle

par.max_radius        =  getParameter(param, {'energy', 'max_radius'}, cutoff_radius );
par.max_radius_factor =  getParameter(param, {'energy', 'radius_factor'}, 10.0);
par.background        =  getParameter(param, {'energy', 'backgournd'}, 0.0);

%inital guess for radii
pari.relative_change    = getParameter(param, {'guess', 'threshold', 'relative_change'}, 0.8); % threshold for relative change in intensity
pari.gradpeak_threshold = getParameter(param, {'guess', 'threshold', 'gradient_peak'},   0.1); % minimal peak height in gradient profile
pari.trough_threshold   = getParameter(param, {'guess', 'threshold', 'trough'},   0.1);        % minimal trough height in itensity profile
pari.background         = getParameter(param, {'guess', 'threshold', 'background', },   0.0);  % definately stop ray below this intensity


% sizes, rays, etc 
[h, w] = size(image);

phi = 0:2*pi/nrays: (2*pi - 2*pi/nrays);
rays = [cos(phi); sin(phi)];
r0 = [h; w] / 2;


% initial guess for radii
for ri = nrays:-1:1
   radius(ri) = guessRadius(image, imagegrad, r0, rays(ri), radius_cutoff, pari);
end


%debug
showState(image, imagegrad, rays, radius);

% rund gradient descent on energy

for t = 1:max_steps
   gradE = gradEnergy(image, imagegrad, rays, pos0, radius, delta_radius, par);
   
   if sum(abs(gradE)) < stop_precision
      break
   else
      radius = radius - time_scale * gradE;
   end
   
   showState(image, imagegrad, rays, radius);
   title(sprintf('time: %g  norm(gradE): %g', t, sum(abs(gradE))));
   
end



%segments = zeros(size(image));

segments = gradE;

end









%% Subfunctions


function guessRadius(image, imagegrad, r0, rays, radius_cutoff, param)

% initialize
nrays = size(rays, 2);

r1 = r0 + radius_cutoff * rays;
x = [r0(1,:)', r1(1,:)'];
y = [r0(2,:)', r1(2,:)'];

nraypixel = 2 * radius_cutoff; % might want to oversample a bit

profiles = zeros(nraypixel, nlines);
for i = 1:nrays
   profiles(:, i) = improfile(imgage, x(i,:), y(i,:), nraypixel);
end

profilesgrad = zeros(nraypixel, nlines);
for i = 1:nrays
   profilesgrad(:, i) = improfile(imagegrad, x(i,:), y(i,:), nraypixel);
end

% relative intensities
relprofiles = profiles ./ repmat(profiles(1,:), nraypixel,1);
relprofilesgrad = profilesgrad ./ repmat(profilesgrad(1,:), nraypixel,1);

%
% find estimate for positions on poly
%

abschange = abs(diff(profiles));
abschangegrad = abs(diff(profilesgrad));

relchange = abs(relprofiles - 1);
relchangegrad = abs(relprofilesgrad - 1);


% bound by backgournd
for i=nrays:-1:1
   rm = find(profiles(:,i)< param.background, 1, 'first');
   if isempty(rm)
      radius_min(i) = nraypixel;
   else
      radius_min(i) = rm;
   end
end





% take a initial guess by change larger than 1000;



for i = nlines:-1:1
   pst1 = find(change(:,i) > changethreshold,1, 'first');
   if isempty(pst1)
      pst(i) = maxlinelength;
   else
      pst(i) = pst1;
   end
end
pos = pst;


% check peaks in gradient
for i = nlines:-1:1
   pks1 = findPeaks(profilesgrad(:,i), gradpeakthreshold, 1);
   if isempty(pks1)
      pks(i) = maxlinelength;
   else
      pks(i) = pks1(1);
   end
end  

% check for troughs in intensity
for i = nlines:-1:1
   trs1 = findTroughs(relprofiles(:,i), troughthreshold, 1);
   if isempty(trs1)
      trs(i) = maxlinelength;
   else
      trs(i) = trs1(1);
   end
end  





for i = 1: nlines
   if pos(i) > pks(i) 
      pos(i) = pks(i);
   end
   if pos(i) > posmin(i)
      pos(i) = posmin(i);
   end
   if pos(i) > trs(i)
      pos(i) = trs(i);
   end
end

pos = vertcat(pos,pos) .* rayvec + r0;


figure(2)
subplot(2,2,1)
hold on
plot(posmin, profiles(sub2ind(size(profiles),posmin,1:nlines)), '*');
plot(trs, profiles(sub2ind(size(profiles),trs,1:nlines)), 'b*');

subplot(2,2,2)
hold on
plot(pks, profilesgrad(sub2ind(size(profilesgrad),pks,1:nlines)), '*');

figure(3)
subplot(2,2,1)
hold on
plot(pst, change(sub2ind(size(change),pst,1:nlines)), '*');


% polygons

if exist('ipoly','var')
   ipoly.delete;
   ipolygrad.delete;
end

ipoly = impoly(ax(1), pos');
ipolygrad = impoly(ax(2), pos');

ipoly.addNewPositionCallback(@ipolygrad.setPosition);
ipolygrad.addNewPositionCallback(@ipoly.setPosition);



end





function gradE = gradEnergy(image, imagegrad, rays, r0, radius, delta_radius, param)

%initialize
nrays = size(rays,2);

pos = [radius; radius] .* rays + repmat(r0,1,nrays);
pos_lo = [pos(:, end) pos(:, 1:end-1)];
pos_hi = [pos(:, 2:end) pos(:, 1)];

pos_new = pos + delta_radius * rays;

diffvec = pos_hi - pos;

diffvec_new_lo = pos_new - pos_lo;
diffvec_new_hi = pos_hi - pos_new;

dist2 = sum(diffvec.^2, 1);
dist = sqrt(dist2);

dist2_new_lo = sum(diffvec_new_lo.^2, 1);
dist2_new_hi = sum(diffvec_new_hi.^2, 1);
dist_new_lo = sqrt(dist2_new_lo);
dist_new_hi = sqrt(dist2_new_hi);

%normdiffvec = 1/dist;
%normdiffvec = [normdiffvec; normdiffvec] .* diffvec;

gradE = zeros(1,nrays);


%visualize
linecl(pos(1,:), pos(2,:))
linecl(pos_new(1,:), pos_new(2,:), 'Color', 'r')


%% accumulate energy terms

%%%
%%% engergy due to deviation form perfect polygon
%%%

if param.weight_angle > 0

   angle = vector_angle([diffvec(:, end) diffvec(:, 1:end-1)], diffvec);
   angle_lo = [angle(end) angle(1:end-1)];
   angle_hi = [angle(2:end) angle(1)];

   angle_delta = vector_angle(diffvec_new_lo, diffvec_new_hi) - angle;
   angle_delta_lo =  vector_angle([diffvec(:, end-1:end) diffvec(:,1:end-2)], diffvec_new_lo) - angle_lo;
   angle_delta_hi =  vector_angle(diffvec_new_hi, [diffvec(:, 2:end) diffvec(:, 1)]) - angle_hi;

   sin_angle = sin(angle - param.target_angle);
   sin_angle_lo = [sin_angle(end) sin_angle(1:end-1)];
   sin_angle_hi = [sin_angle(2:end) sin_angle(1)];

   cos_angle = cos(angle - param.target_angle);
   cos_angle_lo = [cos_angle(end) cos_angle(1:end-1)];
   cos_angle_hi = [cos_angle(2:end) cos_angle(1)];

   % energy: weight_angle * sum( sin(angle - target_angle).^2 )
   gradE = gradE + param.weight_angle * 2 * ( angle_delta    .* sin_angle    .* cos_angle ...  
                                            + angle_delta_lo .* sin_angle_lo .* cos_angle_lo ...
                                            + angle_delta_hi .* sin_angle_hi .* cos_angle_hi);

end
                                
%%%
%%% surface tension
%%%

if param.weight_surface > 0
   
dist2_delta_lo = dist2_new_lo - [dist2(end) dist2(1:end-1)];
dist2_delta_hi = dist2_new_hi - dist2;

surface = sum(dist);
surface_delta = sqrt(sum(dist2) + dist2_delta_lo + dist2_delta_hi) - surface;

% energy: weight_surface * (surface - target_surface)^2
gradE = gradE + param.weight_surface * 2 * surface_delta * (surface - param.target_surface);

end

%%%
%%% radial penalties
%%%

if param.weight_radius > 0
   
% energy: weight_radius * sum((radius - target_radius).^2)
gradE = gradE +  param.weight_radius * 2 * delta_radius * (radius - param.target_radius);

end

if param.weight_max_radius > 0
   
% energy: weight_max_radius * exp(max_radius_factor * (radius - max_radius)

gradE = gradE +  param.weight_max_radius * param.max_radius_factor * exp(param.max_radius_factor * (radius - param.max_radius)) * delta_radius;

end


%%%
%%% area penalty
%%%

if param.weight_area > 0 
   
area = polyarea(pos(1,:), pos(2,:));

for r =narys:-1:1
   pos2 = pos;
   pos2(:, r) = pos_new(:, r);
   area_delta(r)  = polyarea(pos2(1,:), pos2(2,:)) - area;
end

%energy: weight_area * (area - target_area)^2
gradE = gradE + param.weight_area * 2 * area_delta * (area - param.target_area);

end

%%%
%%% gradient coverage
%%%

if param.weight_gradient > 0

   for r=nrays:-1:1

      r_hi = r+1;
      if r_hi == nrays + 1; r_hi = 1; end
      r_lo = r-1;
      if r_lo == 0; r_lo = nrays; end

      np = (round(dist(r))+1) * 2;  % over sample by factor of 2 to make the change smooth
      intgrad_hi(r) = sum(improfile(imagegrad, pos(1, [r r_hi]), pos(2, [r r_hi]), np)) * dist(r) / ( np);

      np = (round(dist_new_lo(r))+1) * 2;   
      intgrad_new_lo(r) = sum(improfile(imagegrad, [pos(1, r_lo), pos_new(1, r)], [pos(2,r_lo) pos_new(2, r)], np)) * dist_new_lo(r) / np;  

      np = (round(dist_new_hi(r))+1) * 2;   
      intgrad_new_hi(r) = sum(improfile(imagegrad, [pos_new(1, r), pos(1, r_hi)], [pos(2,r) pos_new(2, r_hi)], np)) * dist_new_hi(r) / np;  

   end

   intgrad_lo = [intgrad_hi(end) intgrad_hi(1:end-1)];

   intgrad_delta = (intgrad_new_hi - intgrad_hi) + (intgrad_new_lo - intgrad_lo);

   %energy: - weight_gradient * intgrad
   gradE = gradE - param.weight_gradient * intgrad_delta;

end

%%%
%%% intensity variation
%%%

if param.weight_intensity > 0 && param.weight_variation > 0    
   intensities = image(polypixel(pos, h, w));

   for r = nrays:-1:1

      pos2 = pos;
      pos2(:, r) = pos_new(:, r);

      intensities_new{r} = image(polypixel(pos2, h, w));
   end      

   if param.weight_variation > 0
      cv = std(intensities) / mean(intensities);
      cv_new = cellfun(@(x) std(x)/mean(x), intensities_new);

      % energy: weight_variation * CV
      gradE = gradE + weight_variation * (cv_new - cv);
   end

   if param.weight_intensity > 0
      intintensity = sum(intensities);
      intintensity_new = cellfun(@sum , intensities_new);

      % energy: - weight_intensity * intintensity
      gradE = gradE - param.weight_intensity * (intintensity_new - intintensity);
   end
end

end % gradEnergy



% calculates anlge between two succesive vector r0 and r1
function angle = vector_angle(r0, r1)
   dotv = sum(r0 .* r1, 1);
   detv = r0(1,:) .* r1(2,:) - r0(2,:) .* r1(1,:);

   % angle = atan2(detv, dotv) + pi;  % outer angle
   angle = pi - atan2(detv, dotv);    % inner angle               
end

function pixel = polypixel(r, h, w)
   mask = poly2mask(r(1,:), r(2,:), h, w);
   pixel = find(mask);
end




% visualization for debuggin

function showState(image, imagegrad, rays, radius, max_radius)




   figure(52)
   pos =  [radius; radius] .* rays;
   ax(1)= imsubplot(image);
   imshow(image);
   ax(2) = imsubplot(imagegrad);
   imshow(imagegrad);
   linkaxes(ax, 'xy')
   ln = linecl(pos(1,:), pos(2,:));
   
   
   
   % draw rays and profiles 

col = linspace(0,1, nlines+1);
col = col(1:end-1);
col = hsv2rgb([col', ones( nlines,1), ones(nlines,1)]);

% plot profiles along rays

figure(1)
hold on
for s= 1:2
   subplot(1,2,s)
   for i=1:nlines
      line(x(i,:), y(i,:), 'Color', col(i,:))
   end
end
hold off
end




function showProfiles()


figure(2)
clf
ax2(1) = subplot(2,2,1);
set(gcf,'DefaultAxesColorOrder',col)
plot(profiles)

ax2(2) = subplot(2,2,2);
set(gcf,'DefaultAxesColorOrder',col)
plot(profilesgrad)

ax2(3) = subplot(2,2,3);
set(gcf,'DefaultAxesColorOrder',col)
plot(relprofiles)

ax2(4) = subplot(2,2,4);
set(gcf,'DefaultAxesColorOrder',col)
plot(relprofilesgrad)


linkaxes(ax2, 'x');



% plot
figure(3)
clf
ax3(1) = subplot(2,2,1);
set(gcf,'DefaultAxesColorOrder',col)
plot(change)
title('abs change from start')

ax3(2) = subplot(2,2,2);
set(gcf,'DefaultAxesColorOrder',col)
plot(changegrad)
title('abs grad change from start')

ax3(3) = subplot(2,2,3);
set(gcf,'DefaultAxesColorOrder',col)
plot(relchange)
title('abs rel change')

ax3(4) = subplot(2,2,4);
set(gcf,'DefaultAxesColorOrder',col)
plot(relchangegrad)
title('abs rel grad change')

linkaxes(ax3, 'x');












end

