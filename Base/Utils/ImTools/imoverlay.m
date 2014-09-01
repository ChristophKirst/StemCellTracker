function out = imoverlay(in, mask, color, intensity, rescale)
%
% out = imoverlay(in, mask, color, intensity, rescale)
%
% description:
%    replace pixels in image specifed by mask with specified color
%    if intesnity flag is true the intensity in the new color matches the old one
%
% input:
%    in        input image, 2d, 3d, gray or color in form pq, pqc, pql or pqlc
%    mask      mask where pixels should be changed
%    intensity (optional) make masked pixels intensity the one of the old pixels (false)
%    rescale   (optional) rescale original image intensitites (true)
%
% output:
%    out      overlayed image
%

% based on imoverlay by Steven L. Eddins, The MathWorks, Inc.

if nargin < 3 || isempty(color)
    color = [1 0 0];
end
color = imcolorspec2rgb(color);

if nargin < 4
   intensity = false;
end
if nargin < 5
   rescale = false;
end

mask = (mask ~= 0);

dim = ndims(mask);
idim = ndims(in);

% image class unit8 or unit16
if intensity
   if rescale
      in = im2uint16(double(in) / max(in(:))) / 256; 
   else
      in = im2uint16(in) / 256;
   end
   col = im2uint16(color) / 256;
else
   if rescale
      in = im2uint8(double(in) / max(in(:))); 
   else
      in = im2uint8(in); 
   end
   col = im2uint8(color);
end


% Initialize the red, green, and blue output channels.
if dim == 2
   if idim == 2
      % Input is grayscale.  Initialize all output channels the same.
      out_red   = in;
      out_green = in;
      out_blue  = in;
   elseif idim == 3
      % Input is RGB truecolor.
      out_red   = in(:,:,1);
      out_green = in(:,:,2);
      out_blue  = in(:,:,3);
   else
      error('imoverlay: input image has inconsistent dimension');
   end
   
elseif dim ==3
    if idim == 3
      % Input is grayscale.  Initialize all output channels the same.
      out_red   = in;
      out_green = in;
      out_blue  = in;
      
    elseif idim == 4
      % Input is RGB truecolor.
      out_red   = in(:,:,:,1);
      out_green = in(:,:,:,2);
      out_blue  = in(:,:,:,3);
   else
      error('imoverlay: input image has inconsistent dimension');
   end
   
else
   error('imoverlay: input mask has inconsistent dimension');
end

% replace masked pixels
if intensity
   if dim == 2
      if idim == 2
         inm = in(mask);
      else  % idim == 3
         inm = sqrt(sum(double(in(mask)).^2, 3));
      end
   else % dim == 3
      if idim == 3
         inm = in(mask);
      else % idim == 4
         inm = sqrt(sum(double(in(mask)).^2, 4));
      end
   end
   
   out_red(mask)   = inm * col(1);
   out_green(mask) = inm * col(2);
   out_blue(mask)  = inm * col(3);
   
   out_red(~mask)   = out_red(~mask) * 256;
   out_green(~mask) = out_green(~mask) * 256;
   out_blue(~mask)  = out_blue(~mask) * 256;
   
   out_red   = im2uint8(out_red);
   out_green = im2uint8(out_green); 
   out_blue  = im2uint8(out_blue); 
   
else  % intensity == false
   out_red(mask)   = col(1);
   out_green(mask) = col(2);
   out_blue(mask)  = col(3);
end


% combine colors to pqlc format
if dim == 2
   out = cat(3, out_red, out_green, out_blue);
else
   out = cat(4, out_red, out_green, out_blue);
end