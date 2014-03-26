function stats = imstatistics(labeledimage, stats, statnames, image)
%
% stats = imstatistics(labeledimage, statnames, image)
% stats = imstatistics(labeledimage, stats, statnames, image)
%
% description:
%     add standard statistics fields to the struct stats
%     this function is used to avoid calculating statists more than once
%
% input:
%     labeledimage  labeled image, assumed to be relabled with regions 1:max(labeledimage(:))
%     statsname     name or cell of the statistics fields to calculate
%                   regoinporps based statistics:
%                   {'Volume', 'PixelIdxList', 'PixelList', 'PixelValues', 'BoundingBox', 'Centroid', 'WeightedCentroid'}
%                   measurement statistics (additiona image input needed)
%                   {'MaxInstensity', 'MinInstensity', 'MeanIntensity', 'MedianIntensity', 'TotalIntensity', 'StdIntensity', 'VarIntensity'}
%                   if 'XXXIntenisty' function lower(XXX) is applied to the intensity values of each labled region
%                   additional statistics and measures:
%                   {'UltimateErosion', 'SurfacePixelIdxList'} Todo: {'Surface', 'Ellipsoid', 'MaxDiameter', 'MinDiameter'}                 
%     stats         (optional) previously calcualted statistics to avoid recalculation of statistics for existing fields
%     image         (optional) additional intensity image
%
% output:
%     stats         statistics struct with the requested entries
%
% note:
%     imstatistics corrects for regionprops not using pixel coordinates but x,y instead, BoundingBox is calcualted as in imboundingbox
%     'Volume' is used instead of 'Area'
%
% See also: imsurface, impixelsurface, regionprops

%define all statistcs names
regionpropnames =  {'Volume', 'PixelIdxList', 'PixelList', 'PixelValues', 'BoundingBox', 'Centroid', 'WeightedCentroid'};
regionpropnamesimage = {'WeightedCentroid'};
intensitypropnames = {'MaxIntensity', 'MinIntensity', 'MeanIntensity', 'MedianIntensity', '*Intensity'};
specialpropnames =  {'UltimateErosion', 'SurfacePixelIdxList'}; %  'Surface', 'Ellipsoid'
propnames = {regionpropnames{:} regionpropnamesimage{:} intensitypropnames{:} specialpropnames{:}}; %#ok<CCAT>

nlabel = max(labeledimage(:));
dim = ndims(labeledimage);

%get all statistics names
if nargin < 2
   stats = repmat(struct, nlabel, 1);
   statnames = 'default';
end

if isempty(stats)
   stats = repmat(struct, nlabel, 1);
end

if ~isstruct(stats)
   if nargin > 3
      error('imstatistics: inconsistent input')
   end

   if nargin > 2
      image = statnames;
   else
      image  = [];
   end
   statnames = stats;
   stats = repmat(struct, nlabel, 1);
else
   if nargin < 4
      image = [];
   end
end

if ~iscell(statnames)
   statnames = {statnames};
end

if ~iscellstr(statnames)
   error('imstatistics: expect statistic name or cell list of statistics names')
end

if any(ismember(statnames, 'all'))
   statnames = propnames;
end

if any(ismember(statnames, 'default'))
   statnames = {statnames{:}, 'Volume', 'PixelIdxList'}; %#ok<CCAT>
end

statnames(strcmp(statnames, 'Area')) = {'Volume'};

statnames = unique(statnames);
statnames = statnames(~any(cell2mat(cellfun(@(y) strcmp(statnames,y), {'default', 'all'}, 'UniformOutput', 0)'),1));

if length(stats) ~= nlabel
   error('imstatistics: inconsistent object sizes!');
end


%remove existing stats from statnames
fn  = fieldnames(stats);
if ~isempty(fn) 
   statnames = statnames(~any(cell2mat(cellfun(@(y) strcmp(statnames,y), fn, 'UniformOutput', 0)'),1));
end

if isempty(statnames)
   return
end

if any(ismember(statnames, intensitypropnames)) || any(ismember(statnames, regionpropnamesimage))
   if isempty(image)
      error('imstatistics: intensity statistics need intensity image!');
   end
end



%%% regoinprops based statistics

%determine regionprop properties
regprops = {};
if any(ismember(statnames, 'Volume'))
   regprops = [regprops, {'Area'}];
end

idx = ismember(statnames, regionpropnames);
regprops = [regprops, statnames(idx)];

idx = ismember(statnames, regionpropnamesimage);
regprops = [regprops, statnames(idx)];

if any(ismember(statnames, intensitypropnames))
   regprops = [regprops, {'PixelIdxList'}];
end

%if any(ismember(statnames, specialpropnames(4:6)))
%   regprops = [regprops, {'PixelIdxList'}];
%end

%if any(ismember(statnames, specialpropnames(1:3)))
if any(ismember(statnames, specialpropnames(1:2)))
   regprops = [regprops, {'BoundingBox'}];
end


if ~isempty(regprops)
   if isempty(image)
      regstats = regionprops(labeledimage, regprops{:});
   else
      regstats = regionprops(labeledimage, image, regprops{:});
   end

   %Area -> Volume
   if isfield(regstats, 'Area')
      [regstats.('Volume')] = regstats.('Area');
      regstats = rmfield(regstats,'Area');
   end
   
   %Centroid: correct for space to pixel coordinates x,y -> p,q
   if isfield(regstats, 'Centroid')
      dat = {regstats.Centroid};
      if length(regstats(1).Centroid) == 2
         dat = cellfun(@(x) x([2,1])', dat, 'UniformOutput', false);
      else
         dat = cellfun(@(x) x([2,1,3])', dat, 'UniformOutput', false);
      end
      [regstats.Centroid] = dat{:};
   end
   
   %WeightedCentroid: correct for space to pixel coordinates x,y -> p,q
   if isfield(regstats, 'WeightedCentroid')
      dat = {regstats.WeightedCentroid};
      if length(regstats(1).WeightedCentroid) == 2
         dat = cellfun(@(x) x([2,1])', dat, 'UniformOutput', false);
      else
         dat = cellfun(@(x) x([2,1,3])', dat, 'UniformOutput', false);
      end
      [regstats.WeightedCentroid] = dat{:};
   end
   
   
   %BoundingBox: correct to imboundingbox standard   
   if isfield(regstats, 'BoundingBox')
      
      dat = reshape([regstats.BoundingBox], dim * 2, []);
      ee = [2 1 3]; ee = ee(1:dim);
      ee = [ee, ee + dim];
      dat = dat(ee, :);
      
      dat(1:dim, :) = round(dat(1:dim,:) + 0.5); % correct for x,y exchange and shift to full pixel
      ee = (dim+1):(2 * dim);
      dat(ee, :) = dat(1:dim, :) + dat(ee, :) - 1;
    
      dat = num2cell(dat,1);
      
      [regstats.BoundingBox] = dat{:};
   end
   
   % add fields to stats
   for sn = fieldnames(regstats)
      [stats.(sn{1})] = regstats.(sn{1});
   end
   
   
end
   



%%% Intensity measurements

idx = find(~cellfun(@isempty, strfind(statnames, 'Intensity')));
for i = idx
   fn = statnames{i};
   fn = lower(fn(1:end-9));
   fun = eval(['@' fn]);
   
   for l = 1:length(stats)
      stats(l).(statnames{i}) = fun(image(stats(l).PixelIdxList));
   end
end


%%% Specialized measurements

% {'UltimateErosion', 'SurfacePixelIdxList', 'Surface', 'Ellipsoid'};
 
 
if any(ismember(statnames, specialpropnames))

   ue = any(ismember(statnames, 'UltimateErosion'));
   ps = any(ismember(statnames, 'PixelSurface'));
   %su = any(ismember(props, 'Surface')); 
   %el = any(ismember(props, 'Ellipsoid')); 
   
   isize = size(labeledimage);
   dim = length(isize);
   
   for l = 1:nlabel
      % reduce calculation to bounding box
      bmin = stats(l).BoundingBox;
      bmax = bmin((dim+1):(2*dim));
      bmin = bmin(1:dim);
      bmax = bmax';
      bmin = bmin';

      bmin = max(bmin - 1, 1);
      bmax = min(bmax + 1, isize);
      obj = imextract(labeledimage, bmin, bmax);
      obj = (obj == l);
      total(obj)
      
      %UltimateEroison
      if ue
         ulter = imfind(bwulterode(obj));
         size(ulter)
         ulter = ulter + repmat(bmin, size(ulter, 2), 1) - 1;
         stats(l).UltimateErosion = ulter;
      end
      
      %PixelSurface
      if ps
         pixsurf = impixelsurface(obj);
         stats(l).SurfacePixelIdxList = find(pixsurf == l) + bmin - 1;
      end
       
      %Surface 2d ?
%       if su
%          % call isosurface / isonormals
%          [fc,vc] = isosurface(obj, 0.5);
%          nr = isonormals(obj, vc);
% 
%          % correct for x,y exchange and assign outputs
%          vc = vc(:, [2 1 3]);
%          stats(l).Surface = {vc + repmat(bmin, size(v,1), 1) - 1, fc,  nr(:,[2 1 3])};
%       end
         
      
      %Ellipsoid 2d ?
   end
end  % sepcial statistics


end
   
   
