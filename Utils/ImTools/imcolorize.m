function image = imcolorize(image, colmap)
%
% image = imcolorize(image, colmap)
%
% description:
%    colorizes the label image
%
% input:
%    image  matrix of labels
%    colmap color map 
%
% output:
%    image  colorized image
%
% See also:
%    label2rgb

if nargin < 2
   colmap = @jet;
end

if ndims(image) == 2 %#ok<ISMAT>
   image = label2rgb(image, colmap, 'k', 'shuffle');
   image = double(image) / double(max(image(:)));
else
   image = label2rgb3d(image, colmap, 'k', 'shuffle');
   image = double(image) / double(max(image(:)));
end
