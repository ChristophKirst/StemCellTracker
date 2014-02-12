function out = imclip(image, min, max)
%
% out = imclip(image, min, max)
%
% description:
%     clips the intensities in the image to lie withing [min,max]
%
% input:
%    image   the image to be clipped
%    min     minimal intensity
%    max     maximal intensity
%
% output:
%    out     clipped intensity image
%
% See also: imcontrast

if nargin == 1
   min = 0;
   max = 0;
elseif nargin == 2
   if length(min) == 2
      max = min(2);
      min = min(1);
   else
      max = min;
      min = 0;
   end
end

out = image;
out(out>max) = max;
out(out<min) = min;

end