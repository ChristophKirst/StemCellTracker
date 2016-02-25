function [imgpost, stats] = postProcessSegments(imglab, varargin)
%
% [imgpost, stats] = postProcessSegments(imglab, img, param)
%
% input:
%    imglab   labeled image (assume labels are from 1:nalabels, use imrelabel if not)
%    img      (optional) intensity img
%    param    (optional) parameter struct with entries
%             .volume.min            minimal volume to keep (0)
%             .volume.max            maximal volume to keep (Inf)
%             .intensity.min         minimal mean intensity to keep (-Inf)
%             .intensity.max         maximal mean intensity to keep (Inf)
%             .intensity.median.min  minimal median intensity (-Inf)
%             .intensity.median.max  maximal median intensity (Inf)
%             .boundaries            clear objects on x,y boundaries (false)
%             .fillholes             fill holes in each z slice after processing segments (true)
%             .relabel               relabel from 1:nimglabnew (true)
%             .smooth                smooth -- todo e.g. using vtk denoising library / java interface
% 
% output:
%    imgpost  post processed imglab
%    stats    calculated statistics for sequential use
%
% note:
%    clearing of boundary objects here ignores touching z boundaries, use imclearborder for this
%    filling holes is done in each slice only, use imfill(..., 'holes') on full 3d img for this
%
% See also: regionprops, imrelabel, imclearborder, imfill


if nargin < 2
    param = [];
    img = [];
else 
    if isstruct(varargin{1}) || ischar(varargin{1})
       param = parseParameter(varargin{:});
       img = [];
    else
       img = varargin{1};
       param = parseParameter(varargin{2:end});
    end
end

volume_min = getParameter(param, {'volume', 'min'}, 0);
volume_max = getParameter(param, {'volume', 'max'}, Inf);

intensity_min = getParameter(param, {'intensity', 'min'}, -Inf);
intensity_max = getParameter(param, {'intensity', 'max'}, Inf);

intensity_median_min = getParameter(param, 'intensity.median.min', -Inf);
intensity_median_max = getParameter(param, 'intensity.median.max', Inf);

intensity_mean_min = getParameter(param, 'intensity.mean.min', -Inf);
intensity_mean_max = getParameter(param, 'intensity.mean.max', Inf);

boundaries = getParameter(param, {'boundaries'}, false);
fillholes  = getParameter(param, {'fillholes'}, true);
relabel    = getParameter(param, {'relabel'}, true);

imgpost = imglab;


% determine stats to calculate
statnames = {};
if volume_min > 0 || volume_max < Inf
   statnames{end+1} = 'Volume';
   statnames{end+1} = 'PixelIdxList';
end
if intensity_min > -Inf || intensity_max < Inf
   if isempty(img)
      error('postProcessSegments: intensity img expected');
   end
   statnames{end+1} = 'MeanIntensity';
   statnames{end+1} = 'PixelIdxList';
end
if intensity_median_min > -Inf || intensity_median_max < Inf
   if isempty(img)
      error('postProcessSegments: intensity img expected');
   end
   statnames{end+1} = 'MedianIntensity';
   statnames{end+1} = 'PixelIdxList';
end
if intensity_mean_min > -Inf || intensity_mean_max < Inf
   if isempty(img)
      error('postProcessSegments: intensity img expected');
   end
   statnames{end+1} = 'MeanIntensity';
   statnames{end+1} = 'PixelIdxList';
end


if ~isempty(statnames)
    if isempty(img)
        stats = imstatistics(imglab, statnames);
    else
        stats = imstatistics(imglab, statnames, img);
    end
end

if volume_min > 0 || volume_max < Inf
   idx = find([stats.Volume] > volume_max);
   for i = idx
      imgpost(stats(i).PixelIdxList) = 0;  
   end
   idx = find([stats.Volume] < volume_min);
   for i = idx
      imgpost(stats(i).PixelIdxList) = 0;  
   end
end

if intensity_min > -Inf || intensity_max < Inf
   
   %for i = length(stats):-1:1
   %   mi{i} = mean(img(stats(i).PixelIdxList));
   %end
   %[stats.MeanIntensity] = mi{:};
   
   idx = find([stats.MeanIntensity] < intensity_min);
   for i = idx
      imgpost(stats(i).PixelIdxList) = 0;  
   end
   idx = find([stats.MeanIntensity] > intensity_max);
   for i = idx
      imgpost(stats(i).PixelIdxList) = 0;  
   end
end

if intensity_median_min > -Inf || intensity_median_max < Inf
      
   idx = find([stats.MedianIntensity] < intensity_median_min);
   for i = idx
      imgpost(stats(i).PixelIdxList) = 0;  
   end
   idx = find([stats.MedianIntensity] > intensity_median_max);
   for i = idx
      imgpost(stats(i).PixelIdxList) = 0;  
   end
end

if intensity_mean_min > -Inf || intensity_mean_max < Inf
      
   idx = find([stats.MeanIntensity] < intensity_mean_min);
   for i = idx
      imgpost(stats(i).PixelIdxList) = 0;  
   end
   idx = find([stats.MeanIntensity] > intensity_mean_max);
   for i = idx
      imgpost(stats(i).PixelIdxList) = 0;  
   end
end



% fill possible holes in each slice
if fillholes
   for s = 1:size(imgpost,3)
      imgpost(:,:,s) = imfill(imgpost(:,:,s), 'holes');
   end
end

% we dont want to clear objects that border in z
if boundaries
   if ndims(imglab) == 3
      pseg = zeros(size(imgpost) + [0 0 2]);
      pseg(:, :, 2:end-1) = imgpost;
      pseg = imclearborder(pseg);
      imgpost = pseg(:,:,2:end-1);
   else
      imgpost = imclearborder(imgpost);
   end
end

if relabel
   imgpost = imrelabel(imgpost);
   
   % todo: optimize using previous stats results !
   if nargout > 1
     if isempty(img)
        stats = imstatistics(imgpost, statnames);
     else
        stats = imstatistics(imgpost, statnames, img);
     end
   end

end

end


