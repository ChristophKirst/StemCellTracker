function implot(image)
%
% implot(image)
%
% description:
%    implot plots image using pixel coordinates [h,w]
%
% input:
%    image     image to show
%
% See also: implot3d, imshow

imshow(permute(image, [2 1 3]))
axis on
xlabel('w'); ylabel('h'); 