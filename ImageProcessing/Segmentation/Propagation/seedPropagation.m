function [segments, distances] = seedPropagation(image, label, mask, param)
%
% [segments, distances] = segmentByPropagation(label, image, param)
%
% description: 
%    propageates a labeled image and stops propagation if distance measure reaches a cutoff
%    this routine is usefull for propagation of a seeds to prevent oversegmentation
%
% input:
%    image      instensity image to be segmented
%    lables     starting lables
%    mask       (optional) restrict propagation to this mask
%    param      (optional) parameter struct with entries
%               .lambda            interpolation between spatial distances and inensity distance (0 = intensity distance only)
%               .cutoff.difference if the difference measure exceeds this value propgation to the pixel is stopped
%               .averaging.ksize   size of averaging for calculation of the intensity differences (1 = local pixel only)
%               .intensity         intensity references to compare intensity to ([] = average intensity over labels, -1 = calculate average using avaraging.ksize box)
%
% output:
%    segments   labels for segments
%    distances  (optional) distances to initial seeds
%

if (nargin < 3 || isempty(mask))
   mask = ones(size(image));
end
if nargin < 4
   param = [];
end

if isstruct(mask) 
   param = mask;
   mask = ones(size(image));
end
mask = mask > 0;

lambda            = getParameter(param, 'lambda', 0);
cutoff_difference = getParameter(param, 'cutoff.difference', Inf);
averaging_ksize   = getParameter(param, 'averaging.ksize', 1);
intensity_refs    = getParameter(param, 'intensity', []);


radius = round(averaging_ksize(1) - 1) /2;
labels = imlabel(label); 

if isempty(intensity_refs) 
   
   intensity_refs = zeros(max(labels), 1);
   for l = labels
      intensity_refs(l) =  mean(image(label == l));
   end
   
elseif numel(intensity_refs) == 1 && intensity_refs(1) == -1  
   
   intensity_refs = zeros(max(labels));
   for l = labels
      intensity_refs(l) = mean(imfiltervalues(image, label == l, fspecial2('average', averaging_ksize)));
   end
end

if max(labels) > length(intensity_refs)
   error('segmentByPropagation: intensity references do not match label number.');
end

intensity_refs = [0; intensity_refs(:)];

dim = ndims(image);
if dim == 2
   [segments, distances] = seedPropagationMEX(image, label, mask, lambda, cutoff_difference, radius, intensity_refs );
elseif dim ==3
   [segments, distances] = seedPropagation3DMEX(image, label, mask, lambda, cutoff_difference,  radius, intensity_refs );   
else
   error('segmentByPropagation: input image must be 2d or 3d gray scale image')
end

end

   