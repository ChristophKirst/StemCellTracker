function rgb = imgray2color(img, col)
%
% rgb = imgray2color(img)
%
% description:
%  turns 2d or 3d grayscale image into rgb using col as base color
%
% input:
%    img  gray scale image
%
% output:
%    rgb  colored image

dim = ndims(image);
if dim < 2 || dim > 3
   error('imgray2color: expect 2d or 3d grayscale image.')
end

col = imcolor2rgb(col);

rgb=zeros([size(img),3]);

if dim ==2
   rgb(:,:,1)=img * col(1);
   rgb(:,:,2)=img * col(2);
   rgb(:,:,3)=img * col(3);
   %rgb = rgb / 255
else   
   rgb(:,:,:, 1)=img * col(1);
   rgb(:,:,:, 2)=img * col(2);
   rgb(:,:,:, 3)=img * col(3);
   %rgb = rgb / 255
end