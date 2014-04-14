function implot(image, varargin)
%
% implot(image, varargin)
%
% description:
%    implot plots image using pixel coordinates [p,q]
%
% input:
%    image     image to show
%    varargin  options for imshow
%
% See also: implot3d, imshow

image = squeeze(image);
imshow(permute(image, [2 1 3]), varargin{:})
axis on
xlabel('p'); ylabel('q'); 
set(gca,'YDir','normal');

%matlab does not display image if one wants to rotate it
%view(-90,90)
