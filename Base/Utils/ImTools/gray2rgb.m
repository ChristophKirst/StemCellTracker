function rgb = gray2rgb(image)
%
% image = gray2rgb(image)
%
% description:
%  turns 2d or 3d grayscale image into rgb
%

%if ~ismatrix(image)
%   return
%end
dim = ndims(image);
if dim < 2 || dim > 3
   error('gray2rgb: expect 2d or 3d grayscale image.')
end

rgb=zeros([size(image),3]);

if dim ==2
   rgb(:,:,1)=image;
   rgb(:,:,2)=rgb(:,:,1);
   rgb(:,:,3)=rgb(:,:,1);
   %rgb = rgb / 255
else   
   rgb(:,:,:, 1)=image;
   rgb(:,:,:, 2)=rgb(:,:,:, 1);
   rgb(:,:,:, 3)=rgb(:,:,:, 1);
   %rgb = rgb / 255
end