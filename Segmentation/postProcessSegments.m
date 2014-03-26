function [postlabel, stats] = postProcessSegments(label, image, param)
%
% [postlabel, stats] = postProcessSegments(label, image, param)
%
% input:
%    label    labeled image (assume labels 1-nlables, imrelabel)
%    image    (optional) intensity image
%    param    (optional) parameter struct with entries
%             .volume.min        minimal volume to keep (0)
%             .volume.max        maximal volume to keep (Inf)
%             .intensity.min     minimal mean intensity to keep (-Inf)
%             .intensity.max     maximal mean intensity to keep (Inf)
%             .boundaries        clear objects on x,y boundaries (false)
%             .fillholes         fill holes in each z slice after processing segments (true)
%             .relabel           relabel from 1:nlabelnew (true)
%             .smooth            smooth -- todo e.g. using vtk denoising library / java interface
% 
% output:
%    postlabel  post processed label
%    stats      calculated statistics
%
% note:
%    clearing of boundary objects here ignores touching z boundaries, use imclearborder for this
%    filling holes is done in each slice only, use imfill(..., 'holes') on full 3d image for this
%
% Seel also: regionstats, imrelabel, imclearborder, imfill

if nargin < 2
    param = [];
    image = [];
end

if nargin == 2
    if isstruct(image)
        param = image;
        image = [];
    end
end

volume_min = getParameter(param, {'volume', 'min'}, 0);
volume_max = getParameter(param, {'volume', 'max'}, Inf);

intensity_min = getParameter(param, {'intensity', 'min'}, -Inf);
intensity_max = getParameter(param, {'intensity', 'max'}, Inf);

boundaries = getParameter(param, {'boundaries'}, false);
fillholes  = getParameter(param, {'fillholes'}, true);
relabel    = getParameter(param, {'relabel'}, true);

postlabel = label;


% determine stats to calculate
statnames = {};
if volume_min > 0 || volume_max < Inf
   statnames{end+1} = 'Area';
   statnames{end+1} = 'PixelIdxList';
end
if intensity_min > -Inf || intensity_max < Inf
   if isempty(image)
      error('postProcessSegments: intensity image expected');
   end
   statnames{end+1} = 'MeanIntensity';
   statnames{end+1} = 'PixelIdxList';
end

if ~isempty(statnames)
    if isempty(image)
        stats = regionprops(label, statnames{:});
    else
        stats = regionprops(label, image, statnames{:});
    end
end

if volume_min > 0 || volume_max < Inf
   idx = find([stats.Area] > volume_max);
   for i = idx
      postlabel(stats(i).PixelIdxList) = 0;  
   end
   idx = find([stats.Area] < volume_min);
   for i = idx
      postlabel(stats(i).PixelIdxList) = 0;  
   end
end

if intensity_min > -Inf || intensity_max < Inf
   
   %for i = length(stats):-1:1
   %   mi{i} = mean(image(stats(i).PixelIdxList));
   %end
   %[stats.MeanIntensity] = mi{:};
   
   idx = find([stats.MeanIntensity] < intensity_min);
   for i = idx
      postlabel(stats(i).PixelIdxList) = 0;  
   end
   idx = find([stats.MeanIntensity] > intensity_max);
   for i = idx
      postlabel(stats(i).PixelIdxList) = 0;  
   end
end

% fill possible holes in each slice
if fillholes
   for s = 1:size(postlabel,3)
      postlabel(:,:,s) = imfill(postlabel(:,:,s), 'holes');
   end
end

% we dont want to clear objects that border in z
if boundaries
   if ndims(label) == 3
      pseg = zeros(size(postlabel) + [0 0 2]);
      pseg(:, :, 2:end-1) = postlabel;
      pseg = imclearborder(pseg);
      postlabel = pseg(:,:,2:end-1);
   else
      postlabel = imclearborder(postlabel);
   end
end

if relabel
   postlabel = imrelabel(postlabel);
   
   if nargout > 1
     if isempty(image)
        stats = regionstats(label, statnames{:});
     else
        stats = regionstats(label, image, statnames{:});
     end
   end

end

end


