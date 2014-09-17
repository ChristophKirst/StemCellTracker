function rgb = imgray2color(img, colspec)
%
% rgb = imgray2color(img)
% rgb = imgray2color(img, colspec)
% rgb = imgray2color(img, param)
%
% description:
%  turns 2d or 3d grayscale image into rgb using col as base color
%
% input:
%    img     gray scale image
%    colspec (optional) color specification, rgb vector or color name ('r')
%    param   (optional) prameter struct with entries
%            .color.map     use this colormap (colormap)
%            .color.scale   scale color data (true)
%
% output:
%    rgb  colored image
%

dim = ndims(img);
if dim < 2 || dim > 3
   error('imgray2color: expect 2d or 3d grayscale image.')
end

if nargin < 2
   colspec = 'w';
end

if isstruct(colspec)
   colmap = getParameter(colspec, 'color.map', 'default');
   if ischar(colmap)
      colmap = imcolormap(colmap);
   end
   
   colscl = getParameter(colspec, 'color.scale', true);
   si = size(img);
   nm = size(colmap,1);
   
   if colscl
      mx = max(img(:));
      mi = min(img(:));
      if mx > mi
         img = (double(img) -mi)/(mx - mi);
      else
         img = zeros(size(img));
      end      
      rgb = colmap(floor(img(:) * (nm-1)) + 1   , :);

   else
      img(img<1) = 1;
      img(img>nm) = nm;
      rgb = colmap(img(:) , :);
   end
   
   rgb = reshape(rgb, [si, 3]);

else % colspec not struct -> single color

   col = imcolorspec2rgb(colspec);
   rgb = zeros([size(img),3]);
   
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
   
end