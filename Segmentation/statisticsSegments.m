function stats = statisticsSegments(image, label, props)
%
% stats = statisticsSegments(image, label, props)
%
% description:
%      generates statistics for labels in labeled image label using the underlying image data
%
% input:
%      image   image to measure from
%      label   labeled image indicating regions to calc statistics for
%      props   (optional) cell array of names of statistical properties to calculate ([] = {} = 'all')
%
% output:
%      stats   struct array with results
%
% note: many regionprops properties do not work in 3D, first two coordinates are interchanged
%
% See also: regionprops

propnames = {'Volume', 'Centroid', 'BoundingBox', ...
             'MaxIntensity', 'MinIntensity', 'MeanIntensity', 'MedianIntensity', ...
             'UltimateErosion', 'PixelSurface', 'Surface'};

if nargin < 3
   props = 'default';
end

if ~iscell(props)
   props = {props};
end

if any(ismember(props, 'all'))
   props = propnames;
end

if any(ismember(props, 'default'))
   props = propnames(1:7);
end

props = unique(props);

%determine region prop properties
regprops = {};
if any(ismember(props, 'Volume'))
   regprops = [regprops, {'Area'}];
end

for k = 2:3
   if any(ismember(props, propnames{k}))
      regprops = [regprops, propnames(k)]; %#ok<AGROW>
   end
end

if any(ismember(props, propnames(4:7)))
   regprops = [regprops, {'PixelIdxList'}];
end

stats = regionprops(label, regprops{:});

if isempty(stats)
   return
end

% matlab x,y -> h,w
if any(ismember(props, 'Centroid'))
   cent = {stats.Centroid};
   if length(stats(1).Centroid) == 2
      cent = cellfun(@(x) x([2,1])', cent, 'UniformOutput', false);
   else
      cent = cellfun(@(x) x([2,1,3])', cent, 'UniformOutput', false);
   end
   [stats.Centroid] = cent{:};   
end

if any(ismember(props, 'BoundingBox'))
   cent = {stats.BoundingBox};
   if length(stats(1).BoundingBox) == 4
      cent = cellfun(@(x) x([2,1,4,3])', cent, 'UniformOutput', false);
   else
      cent = cellfun(@(x) x([2,1,3,5,4,6])', cent, 'UniformOutput', false);
   end
   [stats.BoundingBox] = cent{:};   
end


%Area -> Volume
if isfield(stats, 'Area')
   [stats.('Volume')] = stats.('Area');
   stats = rmfield(stats,'Area');   
end


%Intensity measurements
idx = find(~cellfun(@isempty, strfind(props, 'Intensity')));
for i = idx
   fn = props{i};
   fn = lower(fn(1:end-9));
   fun = eval(['@' fn]);
   
   for l = 1:length(stats)
      stats(l).(props{i}) = fun(image(stats(l).PixelIdxList));
   end
end


% todo: speed up code below

%label

% if any(ismember(props, propnames(8:10)))
%    %labs = imlabel(label);
%    labs = 1:max(label(:));
%    bb = imlabelboundingboxes(labeledimage); 
%    
%    ue = any(ismember(props, 'UltimateErosion'));
%    ps = any(ismember(props, 'PixelSurface'));
%    su = any(ismember(props, 'Surface')); 
%    
%    for l = labs
%       % reduce calculation to bounding box
%       bmin = bb(l, 1:3);
%       bmax = bb(l, 4:6);
% 
%       bmin = max(bmin - 1, 1);
%       bmax = min(bmax + 1, isize);
%       obj = imextract(labeledimage, bmin, bmax);
%       obj = (obj == l);
%       
%       %UltimateEroison
%       if ue
%          stats(l).UltimateErosion = find(bwulterode(obj)) + bmin - 1;
%       end
%       
%       %PixelSurface
%       if ps
% 
%    for l = labs
%       stats(i).PixelSurface = find(pixsurf == l);
%       i = 1 + 1;
%    end
%       end
% 
%  
% 
%          if ps
%       pixsurf = impixelsurface(label);
%         
%    
%    
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
%    % call isosurface / isonormals
%    [f,v] = isosurface(obj, 0.5);
%    if nargout == 3
%       n = isonormals(obj, v);
%    end
%    
%    % correct for x,y exchange and assign outputs
%    v = v(:, [2 1 3]);
%    vertices{l} = v + repmat(bmin, size(v,1), 1) - 1;
%    faces{l} = f;
%    
%    if nargout == 3
%       normals{l} = n(:,[2 1 3]);
%    end
% 


end