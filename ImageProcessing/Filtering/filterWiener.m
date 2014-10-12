function [img, noise] = filterWiener(img, ksize)
%
% img = filterWiener(img, ksize)
%
% description:
%    wiener filter (for 2d images only)
%
% input:
%    img          image to be filtered
%    ksize        h x w (xl) kernel size
%
% output:
%    img          filtered image

if ismatrix(imgs)
   [img, noise] = wiener2( img, ksize );
else
   error('filterWiener: not implemented for 3d yet');
end

end
   