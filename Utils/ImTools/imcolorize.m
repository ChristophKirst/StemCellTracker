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

image = label2rgb(image, colmap, 'k', 'shuffle');
