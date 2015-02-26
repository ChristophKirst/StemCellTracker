function [segments, distances] = segmentByPropagation(image, label, mask, lambda, ksize, intensity_refs)
%
% [segments, distances] = segmentByPropagation(label, image, mask, lambda)
%
% description: 
%    segments image using propagation algorithm to neighbouring pixels
%    neighbours are included if distance is small enough
%    distance is calucalted as relative intensity change to make algorithm
%    robust against differently illuminated objects
%
% input:
%    image           instensity image to be segmented
%    lables          starting lables
%    mask            restrict propagation to this mask
%    lambda          weight for spatial distances (weight for intensity distances = 1)
%    ksize           width/height of box to calculate intensity differences for (1 = center pixel only)
%    intensity_refs  intensity references ([] = use average over seeds labels, -1 use average over box of size ksize)
%
% output:
%    segments   labels for segments
%    distances  distances to initial seeds
%


if nargin < 3
   mask = ones(size(image));
end
mask = mask > 0;

if nargin < 4
   lambda = 1.0;
end

if nargin < 5
   ksize = 3;
end
radius = round(ksize(1) -1) /2;

if nargin < 6
   intensity_refs = [];
else
   intensity_refs = intensity_refs(:);
end


labels = imlabel(label); 

if isempty(intensity_refs) 
   
   intensity_refs = zeros(max(labels), 1);
   for l = labels
      intensity_refs(l) =  mean(image(label == l));
   end
   
elseif numel(intensity_refs) == 1 && intensity_refs(1) == -1  
   
   intensity_refs = zeros(max(labels));
   for l = labels
      intensity_refs(l) = mean(imfiltervalues(image, label == l, fspecial2('average', kisze)));
   end
end

if max(labels) > length(intensity_refs)
   error('segmentByPropagation: intensity references do not match label number.');
end

intensity_refs = [0; intensity_refs(:)];

%figure(98)
%plot(intensity_refs)

dim = ndims(image);
if dim == 2
   [segments, distances] = segmentByPropagationMEX(image, label, mask, lambda, radius, intensity_refs );
elseif dim ==3
   [segments, distances] = segmentByPropagation3DMEX(image, label, mask, lambda, radius, intensity_refs );   
else
   error('segmentByPropagation: input image must be 2d or 3d gray scale image')
end

end

   