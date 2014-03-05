function objects = label2objects(image, label, time, rescale)
%
% objects = label2objects(label)
%
% description:
%    converts a labeld image to array of Object classes for tracking
%
% input:
%    image   original grayscale image to measure data from (2D/3D)
%    label   labeled image (2D/3D)
%    time    (otional) time point (0)
%    rescale (optional) rescale the coordinates (1)
%    
%
% output:
%    objects   array of Objects each representing one of the labels
%
% See also: Object

if nargin < 2
   time = 0;
end

if nargin < 3
   rescale = 1;
end
rescale = rescale(:)';


prop = regionprops(label, 'Area', 'Centroid', 'PixelIdxList');

%correct for pixel versus space coordinates
cent = {prop.Centroid};
dim = ndims(label);

if length(rescale) < dim
   rescale = rescale * ones(1, ndim);
else
   rescale = rescale(1:dim);
end
   
if dim == 2
   cent = cellfun(@(x) x([2,1]), cent, 'UniformOutput', false);
else
   cent = cellfun(@(x) x([2,1,3]), cent, 'UniformOutput', false);
end

[prop.Centroid] = cent{:};


nobj = length(prop);
for i = nobj:-1:1
   objects(i) = Object(label(prop(i).PixelIdxList(1)), time, rescale' .* prop(i).Centroid', prop(i).Area, median(image(prop(i).PixelIdxList)));
end


