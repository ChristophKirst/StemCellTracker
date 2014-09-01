function out = imclip(image, cmin, cmax)
%
% out = imclip(image, min, max)
%
% description:
%     clips the intensities in the image to lie within [min,max]
%
% input:
%    image   the image to be clipped
%    cmin    minimal intensity ([] = min(image(:)))
%    cmax    maximal intensity ([] = max(image(:)))
%
% output:
%    out     clipped intensity image
%
% See also: imcontrast

if nargin == 1
   cmin = 0;
   cmax = [];
elseif nargin == 2
   if length(cmin) == 2
      cmax = cmin(2);
      cmin = cmin(1);
   else
      cmax = [];
   end
end

if isempty(cmin)
   cmin = min(image(:));
end
if isempty(cmax)
   cmax = max(image(:));
end

out = image;
out(out>cmax) = cmax;
out(out<cmin) = cmin;

end