function image = gray2rgb(image)
%
% image = gray2rgb(image)
%
% description:
%  turns grayscale image into rgb
%

if ~ismatrix(image)
   return
end

[n, m]=size(image);
rgb=zeros(n,m,3);
rgb(:,:,1)=image;
rgb(:,:,2)=rgb(:,:,1);
rgb(:,:,3)=rgb(:,:,1);
image=rgb/255;

end